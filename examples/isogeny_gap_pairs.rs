//! # Isogenous pairs with a large conductor gap
//!
//! The constructor behind
//! `research/isogeny_conductor_gap_ic_20261010/PROTOCOL.md` part B.  For
//! each family it builds two `F_q`-isogenous curves whose endomorphism
//! rings sit at the two ends of an `ℓ`-volcano — the **crater**
//! (`End = O_K`) and the **floor** (`End = Z[π]`) — joined by an explicit
//! chain of `ℓ`-isogenies of total degree `ℓ^h`, the conductor gap, and
//! writes every curve as an explicit `ecbench` workload.
//!
//! * **prime**: CM construction.  `D_K` of class number one, `f_π = ℓ^h`,
//!   and the trace chosen with `(t − f D_K)/2 ≡ 1 (mod ℓ)` so that the
//!   full `ℓ`-torsion is rational on every curve above the floor
//!   (`E[ℓ^k] ⊂ E(F_p)` iff `(π − 1)/ℓ^k ∈ End(E)`).  Vélu descends one
//!   level per step; a step is descending exactly when its codomain's
//!   `j` is new, since distinct levels have distinct endomorphism rings
//!   and hence distinct `j`.  The floor is certified by a cyclic
//!   `ℓ`-part of `E(F_p)`.
//! * **koblitz**: `K_a / F_{2^n}` is the crater of its class
//!   (`End = Z[τ] = O_K`, `D_K = −7`).  For every `ℓ ∣ f_π` the
//!   `ℓ`-torsion is built in `F_{q^m}`, `m = ord_ℓ(t/2)`, with
//!   `binary_torsion_walk`, one descending kernel is taken per level, and
//!   the floor is certified by rank-1 `E[ℓ]` over the same extension.
//! * **binary**: random binary curves whose Frobenius conductor has a
//!   prime factor `ℓ ≥ min_ell`; the vertex found is placed by the rank
//!   of `E[ℓ]` over `F_{q^m}` and the other end of the height-one
//!   volcano is reached by one edge.
//!
//! Every edge is checked: kernel points of order `ℓ` (prime) or a
//! rational kernel trace (binary), `#E′ = #E`, the image of the
//! generator has order `r`, and a planted scalar transports.
//!
//! ```bash
//! cargo run --release --example isogeny_gap_pairs -- --family prime --out DIR
//! cargo run --release --example isogeny_gap_pairs -- --family koblitz --out DIR
//! cargo run --release --example isogeny_gap_pairs -- --family binary --out DIR
//! ```

use std::collections::BTreeSet;
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::Zero;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

use crypto_lib::binary_ecc::IrreduciblePoly;
use crypto_lib::cryptanalysis::binary_torsion_walk::{
    kernel_generators, kernel_of, mult_order, order_over_extension, scalar_eigenvalue,
    torsion_basis, velu_x_from_abscissae, Ext, ExtCurve, Pt, Torsion,
};
use crypto_lib::cryptanalysis::binary_velu::Curve as BinCurve;
use crypto_lib::cryptanalysis::ic_boundary::{
    binary_point_count, koblitz_instance, random_binary_instance, ArtinSchreier, BinaryInstance,
    CountedGroup, GroupOps, PrimeCurve, PrimePoint,
};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use crypto_lib::isogeny::cm::{
    class_number_one_j, cm_discriminant, fundamental_discriminant_and_conductor,
};
use crypto_lib::isogeny::velu::{velu_isogeny_from_kernel, VeluIsogeny};
use crypto_lib::isogeny::volcano::j_invariant;
use crypto_lib::isogeny::SmallCurve;

// ── word arithmetic ─────────────────────────────────────────────────

fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
fn addmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 + b as u128) % p as u128) as u64
}
fn submod(a: u64, b: u64, p: u64) -> u64 {
    addmod(a, p - b % p, p)
}
fn powmod(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1 % p;
    a %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, p);
        }
        a = mulmod(a, a, p);
        e >>= 1;
    }
    r
}
fn invmod(a: u64, p: u64) -> u64 {
    powmod(a, p - 2, p)
}
fn modp_i(v: i128, p: u64) -> u64 {
    let m = v.rem_euclid(p as i128);
    m as u64
}

/// Trial division up to `2^20`, then a primality test on what is left.
fn factor_u64(mut v: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut d = 2u64;
    while d * d <= v && d < (1 << 20) {
        let mut e = 0;
        while v % d == 0 {
            v /= d;
            e += 1;
        }
        if e > 0 {
            out.push((d, e));
        }
        d += if d == 2 { 1 } else { 2 };
    }
    if v > 1 {
        if is_prime_u64(v) {
            out.push((v, 1));
        } else {
            // Composite above the trial bound: record it as is, flagged
            // by exponent 0 so a caller never mistakes it for a prime.
            out.push((v, 0));
        }
    }
    out
}

fn valuation(mut v: u64, ell: u64) -> u32 {
    let mut k = 0;
    while v > 0 && v % ell == 0 {
        v /= ell;
        k += 1;
    }
    k
}

// ── prime fields: the CM ladder ─────────────────────────────────────

struct PrimeLevel {
    curve: PrimeCurve,
    level: u32,
    j: u64,
    generator: PrimePoint,
    aut: u32,
    /// `(x, y)` of the kernel generator of the edge that left this level.
    edge_kernel: Option<(u64, u64)>,
    /// Conductor of `End` at this level, from the chain.
    conductor: u64,
}

struct PrimeCtx {
    p: u64,
    n: u64,
    r: u64,
    ell: u64,
    rng: StdRng,
}

fn small(c: &PrimeCurve) -> SmallCurve {
    SmallCurve {
        name: "cm",
        p: c.p,
        a: c.a,
        b: c.b,
    }
}

fn random_point(c: &PrimeCurve, rng: &mut StdRng) -> PrimePoint {
    loop {
        let x = rng.gen_range(0..c.p);
        let rhs = c.rhs(x);
        if let Some(y) = c.sqrt(rhs) {
            if mulmod(y, y, c.p) == rhs {
                return PrimePoint::affine(x, y);
            }
        }
    }
}

fn pt_mul(c: &PrimeCurve, p: PrimePoint, k: u64) -> PrimePoint {
    let mut ops = GroupOps::default();
    c.mul(&mut ops, p, k)
}
fn pt_add(c: &PrimeCurve, p: PrimePoint, q: PrimePoint) -> PrimePoint {
    let mut ops = GroupOps::default();
    c.add(&mut ops, p, q)
}

/// Does `[n]P = O` hold on `samples` random points?
fn order_holds(c: &PrimeCurve, n: u64, samples: usize, rng: &mut StdRng) -> bool {
    (0..samples).all(|_| pt_mul(c, random_point(c, rng), n).infinity)
}

/// The `ℓ`-part of `E(F_p)`: points of order `ℓ` found, and whether two
/// of them are independent (rank 2: above the `ℓ`-floor) or every one lies
/// on a single line (rank 1: the `ℓ`-floor).
struct EllPart {
    order_ell: Vec<PrimePoint>,
    rank2: bool,
}

fn ell_part_of(c: &PrimeCurve, n: u64, ell: u64, samples: usize, rng: &mut StdRng) -> EllPart {
    let v = valuation(n, ell);
    let m = n / ell.pow(v);
    let mut order_ell = Vec::new();
    for _ in 0..samples {
        let mut q = pt_mul(c, random_point(c, rng), m);
        if q.infinity {
            continue;
        }
        loop {
            let lq = pt_mul(c, q, ell);
            if lq.infinity {
                break;
            }
            q = lq;
        }
        order_ell.push(q);
    }
    let mut rank2 = false;
    if let Some(&q1) = order_ell.first() {
        let line: Vec<PrimePoint> = (0..ell).map(|k| pt_mul(c, q1, k)).collect();
        rank2 = order_ell.iter().any(|q| !line.contains(q));
    }
    EllPart { order_ell, rank2 }
}

fn ell_part(ctx: &mut PrimeCtx, c: &PrimeCurve, samples: usize) -> EllPart {
    ell_part_of(c, ctx.n, ctx.ell, samples, &mut ctx.rng)
}

/// Two independent points of order `ℓ`, or `None`.
fn ell_basis(ctx: &mut PrimeCtx, c: &PrimeCurve) -> Option<(PrimePoint, PrimePoint)> {
    let part = ell_part(ctx, c, 24 + 2 * ctx.ell as usize);
    let q1 = *part.order_ell.first()?;
    let line: Vec<PrimePoint> = (0..ctx.ell).map(|k| pt_mul(c, q1, k)).collect();
    let q2 = *part.order_ell.iter().find(|q| !line.contains(q))?;
    Some((q1, q2))
}

fn velu_from_generator(c: &PrimeCurve, g: PrimePoint, ell: u64) -> VeluIsogeny {
    let pts: Vec<(u64, u64)> = (1..ell)
        .map(|k| {
            let q = pt_mul(c, g, k);
            (q.x, q.y)
        })
        .collect();
    velu_isogeny_from_kernel(&small(c), &pts)
}

fn aut_order(j: u64) -> u32 {
    if j == 0 {
        6
    } else if j == 1728 {
        4
    } else {
        2
    }
}

/// Search `s` for `p = (t² − f²D)/4` prime, `t = fD + 2 + 2ℓs`,
/// `N = ℓ² · c · r` with `r` prime and `c ≤ 4` (`c = 1` where the
/// arithmetic allows it; for `D = −7` the quotient `N/ℓ²` is always even).
/// Returns `(p, t, N, r)`.  The whole window `2^{bits−1} ≤ p < 2^{bits+1}`
/// is scanned from a seed-chosen offset.
fn cm_search(dk: i64, ell: u64, h: u32, bits: u32, seed: u64) -> Option<(u64, i64, u64, u64)> {
    let d = dk as i128;
    // When D ≡ 1 (mod 8) the norm (t² − f²D)/4 with t, f odd is even, so
    // the conductor takes an extra factor 2 and s is kept even so that the
    // full 2-torsion is rational above the 2-floor.
    let two_step = dk.rem_euclid(8) == 1;
    let f = (ell as i128).pow(h) * if two_step { 2 } else { 1 };
    let fd = f * f * (-d); // f²|D|
    let p_lo = 1i128 << (bits - 1);
    let p_hi = 1i128 << (bits + 3);
    let isqrt = |v: i128| -> i128 {
        if v <= 0 {
            return 0;
        }
        let mut x = (v as f64).sqrt() as i128;
        while x * x > v {
            x -= 1;
        }
        while (x + 1) * (x + 1) <= v {
            x += 1;
        }
        x
    };
    let t_lo = isqrt(4 * p_lo - fd).max(1);
    let t_hi = isqrt(4 * p_hi - fd);
    let s_lo = ((t_lo - f * d - 2) / (2 * ell as i128)).max(1);
    let s_hi = (t_hi - f * d - 2) / (2 * ell as i128) + 1;
    if s_hi <= s_lo {
        return None;
    }
    let span = s_hi - s_lo;
    let mut rng = StdRng::seed_from_u64(
        seed ^ (dk as u64) ^ ell.rotate_left(8) ^ u64::from(h) ^ u64::from(bits) << 32,
    );
    let offset = rng.gen_range(0..span);
    for i in 0..span {
        let s = s_lo + (offset + i) % span;
        if s % ell as i128 == 0 || (two_step && s % 2 != 0) {
            continue;
        }
        let t = f * d + 2 + 2 * ell as i128 * s;
        let num = t * t - f * f * d;
        if num % 4 != 0 {
            continue;
        }
        let p = num / 4;
        if p < p_lo || p >= p_hi || p >= (1i128 << 62) {
            continue;
        }
        let p = p as u64;
        if !is_prime_u64(p) || p == ell || p < 5 {
            continue;
        }
        let n = (p as i128 + 1 - t) as u64;
        if valuation(n, ell) != 2 {
            continue;
        }
        let mut r = n / (ell * ell);
        let mut c = 1u64;
        while r % 2 == 0 && c < 32 {
            r /= 2;
            c *= 2;
        }
        let mut c3 = 0;
        while r % 3 == 0 && c3 < 2 {
            r /= 3;
            c3 += 1;
        }
        if r <= ell || r == p || !is_prime_u64(r) || r < (1 << 10) {
            continue;
        }
        return Some((p, t as i64, n, r));
    }
    None
}

/// A curve with `j` and group order `n` over `F_p`, found among the
/// twists by `[n]P = O` on random points.
fn curve_with_j_and_order(p: u64, j: u64, n: u64, rng: &mut StdRng) -> Option<PrimeCurve> {
    let mut cands: Vec<PrimeCurve> = Vec::new();
    if j == 0 {
        for b in 1..64u64 {
            cands.push(PrimeCurve { p, a: 0, b });
        }
    } else if j == 1728 % p {
        for a in 1..64u64 {
            cands.push(PrimeCurve { p, a, b: 0 });
        }
    } else {
        let k = mulmod(j, invmod(submod(1728 % p, j, p), p), p);
        let a = mulmod(3, k, p);
        let b = mulmod(2, k, p);
        cands.push(PrimeCurve { p, a, b });
        // the quadratic twist by every small d
        for d in 2..64u64 {
            cands.push(PrimeCurve {
                p,
                a: mulmod(a, mulmod(d, d, p), p),
                b: mulmod(b, mulmod(mulmod(d, d, p), d, p), p),
            });
        }
    }
    for c in cands {
        let disc = addmod(
            mulmod(4, mulmod(mulmod(c.a, c.a, p), c.a, p), p),
            mulmod(27, mulmod(c.b, c.b, p), p),
            p,
        );
        if disc == 0 {
            continue;
        }
        if j_invariant(&small(&c)) != j {
            continue;
        }
        if order_holds(&c, n, 4, rng) {
            return Some(c);
        }
    }
    None
}

fn prime_explicit_json(c: &PrimeCurve, n: u64, r: u64, g: PrimePoint) -> Value {
    json!({"kind": "prime_explicit", "p": c.p, "a": c.a, "b": c.b, "group_order": n, "r": r, "gx": g.x, "gy": g.y})
}

fn prime_ladder(dk: i64, ell: u64, h: u32, bits: u32, seed: u64) -> Result<Value, String> {
    let started = Instant::now();
    let (p, t, n, r) = cm_search(dk, ell, h, bits, seed).ok_or("no (p, t) found")?;
    let two_step = dk.rem_euclid(8) == 1;
    let steps: Vec<(u64, u32)> = if two_step {
        vec![(ell, h), (2, 1)]
    } else {
        vec![(ell, h)]
    };
    let gap_total: u64 = steps.iter().map(|&(l, e)| l.pow(e)).product();
    let h_total: u32 = steps.iter().map(|&(_, e)| e).sum();
    let mut ctx = PrimeCtx {
        p,
        n,
        r,
        ell,
        rng: StdRng::seed_from_u64(seed ^ p),
    };
    let j_int = class_number_one_j(dk).ok_or("class number one only")?;
    let j0 = modp_i(j_int as i128, p);
    let e0 =
        curve_with_j_and_order(p, j0, n, &mut ctx.rng).ok_or("no crater model with the order")?;
    // Independent check of the trace by the library's point count.
    let cm = cm_discriminant(&small(&e0));
    if cm.trace != t {
        return Err(format!("trace mismatch: built {t}, counted {}", cm.trace));
    }
    let (fund, cond) = fundamental_discriminant_and_conductor(t * t - 4 * p as i64);
    if fund != dk || cond != gap_total as i64 {
        return Err(format!(
            "discriminant mismatch: {fund}·{cond}² against gap {gap_total}"
        ));
    }
    // generator of order r
    let g0 = loop {
        let q = pt_mul(&e0, random_point(&e0, &mut ctx.rng), n / r);
        if !q.infinity && pt_mul(&e0, q, r).infinity {
            break q;
        }
    };
    let planted = 1 + ctx.rng.gen_range(0..r - 1);
    let mut q_planted = pt_mul(&e0, g0, planted);

    let mut levels: Vec<PrimeLevel> = vec![PrimeLevel {
        curve: e0,
        level: 0,
        j: j0,
        generator: g0,
        aut: aut_order(j0),
        edge_kernel: None,
        conductor: 1,
    }];
    let mut visited: BTreeSet<u64> = BTreeSet::new();
    visited.insert(j0);
    let mut edges = Vec::new();
    let mut k = 0u32;
    let mut conductor = 1u64;
    for &(ell, h_ell) in &steps {
        ctx.ell = ell;
        for step in 0..h_ell {
            let steps_left_for_this_ell = h_ell - step;
            let cur = levels.last().unwrap().curve.clone();
            let part = ell_part(&mut ctx, &cur, 24 + 2 * ell as usize);
            if !part.rank2 {
                return Err(format!("level {k}: E[{ell}] not of rank 2"));
            }
            let (q1, q2) = ell_basis(&mut ctx, &cur).ok_or("no ℓ-basis")?;
            let mut gens = vec![q1];
            let mut acc = q2;
            for _ in 0..ell {
                gens.push(acc);
                acc = pt_add(&cur, acc, q1);
            }
            // A kernel descends when its codomain's j is new and the
            // codomain's ℓ-rank is what the next level must have: 2
            // while levels remain, 1 at the ℓ-floor.  The rank test is
            // what rules out a horizontal edge at a non-maximal crater
            // (class number above one), where a new j is not enough.
            let mut chosen: Option<(VeluIsogeny, PrimePoint, u32)> = None;
            let mut skipped_j = Vec::new();
            let mut skipped_rank = Vec::new();
            for (idx, g) in gens.iter().enumerate() {
                if pt_mul(&cur, *g, ell) != PrimePoint::INFINITY || g.infinity {
                    return Err("kernel generator is not of order ℓ".into());
                }
                let iso = velu_from_generator(&cur, *g, ell);
                let jp = j_invariant(&iso.codomain);
                if visited.contains(&jp) {
                    skipped_j.push(jp);
                    continue;
                }
                let cod = PrimeCurve {
                    p,
                    a: iso.codomain.a,
                    b: iso.codomain.b,
                };
                let want_rank2 = steps_left_for_this_ell > 1;
                let cod_part = ell_part_of(&cod, n, ell, 24 + 2 * ell as usize, &mut ctx.rng);
                if cod_part.rank2 != want_rank2 {
                    skipped_rank.push(jp);
                    continue;
                }
                chosen = Some((iso, *g, idx as u32));
                break;
            }
            let (iso, g, idx) = chosen.ok_or(format!("level {k}: no descending kernel"))?;
            if iso.degree != ell {
                return Err(format!("Vélu degree {} ≠ ℓ", iso.degree));
            }
            let next = PrimeCurve {
                p,
                a: iso.codomain.a,
                b: iso.codomain.b,
            };
            if !order_holds(&next, n, 4, &mut ctx.rng) {
                return Err(format!("level {}: #E′ ≠ #E", k + 1));
            }
            let gl = levels.last().unwrap().generator;
            let (gx, gy) = iso.evaluate(gl.x, gl.y).ok_or("generator in the kernel")?;
            let g_next = PrimePoint::affine(gx, gy);
            if !next.is_on_curve(g_next) || !pt_mul(&next, g_next, r).infinity || g_next.infinity {
                return Err(format!(
                    "level {}: image of the generator is not of order r",
                    k + 1
                ));
            }
            let (qx, qy) = iso
                .evaluate(q_planted.x, q_planted.y)
                .ok_or("target in the kernel")?;
            let q_next = PrimePoint::affine(qx, qy);
            if pt_mul(&next, g_next, planted) != q_next {
                return Err(format!("level {}: planted scalar did not transport", k + 1));
            }
            q_planted = q_next;
            let jp = j_invariant(&iso.codomain);
            visited.insert(jp);
            conductor *= ell;
            edges.push(json!({
                "from_level": k, "to_level": k + 1, "ell": ell, "kernel_generator": [g.x, g.y],
                "kernel_index": idx, "codomain": {"a": next.a, "b": next.b}, "codomain_j": jp,
                "skipped_visited_j": skipped_j, "skipped_wrong_rank_j": skipped_rank,
            }));
            levels.last_mut().unwrap().edge_kernel = Some((g.x, g.y));
            levels.push(PrimeLevel {
                curve: next,
                level: k + 1,
                j: jp,
                generator: g_next,
                aut: aut_order(jp),
                edge_kernel: None,
                conductor,
            });
            k += 1;
        }
    }
    // Floor certificate: rank-1 ℓ-part for every ℓ of the chain.
    let floor = levels.last().unwrap().curve.clone();
    let mut floor_ranks = Vec::new();
    for &(l, _) in &steps {
        let part = ell_part_of(&floor, n, l, 24 + 2 * l as usize, &mut ctx.rng);
        if part.rank2 || part.order_ell.is_empty() {
            return Err(format!("floor: E[{l}] is not of rank 1"));
        }
        floor_ranks.push(json!({"ell": l, "rank": 1, "samples": part.order_ell.len()}));
    }
    let floor_cm = cm_discriminant(&small(&floor));
    // Isomorphic model of the floor: (u⁴a, u⁶b), the I-4 control.
    let u = 2 + ctx.rng.gen_range(0..p - 3);
    let u2 = mulmod(u, u, p);
    let u4 = mulmod(u2, u2, p);
    let u6 = mulmod(u4, u2, p);
    let iso_model = PrimeCurve {
        p,
        a: mulmod(floor.a, u4, p),
        b: mulmod(floor.b, u6, p),
    };
    let gf = levels.last().unwrap().generator;
    let g_iso = PrimePoint::affine(mulmod(gf.x, u2, p), mulmod(gf.y, mulmod(u2, u, p), p));
    if !iso_model.is_on_curve(g_iso) || !pt_mul(&iso_model, g_iso, r).infinity {
        return Err("isomorphic model failed".into());
    }
    let curves: Vec<Value> = levels
        .iter()
        .map(|l| {
            json!({
                "role": if l.level == 0 { "crater" } else if l.level == h_total { "floor" } else { "level" },
                "level": l.level,
                "conductor": l.conductor,
                "j": l.j, "aut": l.aut,
                "spec": prime_explicit_json(&l.curve, n, r, l.generator),
            })
        })
        .collect();
    Ok(json!({
        "family": "prime", "d_k": dk, "ell": ell, "h": h, "steps": steps, "gap": gap_total, "bits": bits,
        "p": p, "trace": t, "group_order": n, "r": r, "cofactor": n / r,
        "frobenius_conductor": cond, "fundamental_discriminant": fund,
        "crater": {"j": j0, "aut": aut_order(j0), "end_conductor": 1, "end_evidence": "j = j(O_K), class number one"},
        "floor": {"j": levels.last().unwrap().j, "end_conductor": gap_total,
                   "end_evidence": "rank-1 E[ℓ] for every ℓ | f_π: (π−1)/ℓ ∉ End(E), so End = Z[π]",
                   "ranks": floor_ranks, "cm": floor_cm},
        "planted_scalar": planted,
        "edges": edges,
        "curves": curves,
        "floor_isomorphic_model": {"u": u, "spec": prime_explicit_json(&iso_model, n, r, g_iso)},
        "seconds": started.elapsed().as_secs_f64(),
    }))
}

// ── binary fields ───────────────────────────────────────────────────

fn modulus_mask(irr: &IrreduciblePoly) -> u64 {
    irr.low_terms
        .iter()
        .fold(1u64 << irr.degree, |m, &k| m | (1u64 << k))
}

fn binary_explicit_json(
    n: u32,
    irr: &IrreduciblePoly,
    a: u64,
    b: u64,
    order: u64,
    r: u64,
    g: (u64, u64),
) -> Value {
    json!({"kind": "binary_explicit", "n": n, "modulus": modulus_mask(irr), "a": a, "b": b,
           "group_order": order, "r": r, "gx": g.0, "gy": g.1})
}

/// A point of order `r` above `x` on the curve, or `None`.
fn lift_order_r(curve: &BinCurve, x: u64, r: u64) -> Option<(u64, u64)> {
    let g = curve.group();
    let mut ops = GroupOps::default();
    for pt in curve.points_with_x(x) {
        if g.mul(&mut ops, pt, r).infinity && !pt.infinity {
            return Some((pt.x, pt.y));
        }
    }
    None
}

fn count_order(irr: &IrreduciblePoly, a: u64, b: u64) -> u64 {
    let gf = crypto_lib::cryptanalysis::semaev_decomp::Gf2::new(irr);
    let ash = ArtinSchreier::new(&gf);
    binary_point_count(&gf, &ash, a, b)
}

struct BinaryEdge {
    a6_from: u64,
    a6_to: u64,
    ell: u64,
    m: usize,
    x_image: u64,
    kernels_tried: usize,
    loops_skipped: usize,
}

/// One descending `ℓ`-edge from a rank-2 vertex `(a2, a6)`, carrying
/// abscissa `x_carry` (of a point of order `r`).  `avoid` lists the
/// codomain `a6` values that would be a loop or a backtrack.
fn descend_binary(
    inst: &BinaryInstance,
    ext: &Ext,
    ell: u64,
    order_m: &BigUint,
    a6: u64,
    x_carry: u64,
    avoid: &[u64],
    want_floor: Option<bool>,
    seed: u64,
) -> Result<(BinaryEdge, bool), String> {
    let curve = ExtCurve::new(ext, inst.a, a6);
    let mut rng = StdRng::seed_from_u64(seed ^ a6 ^ ell.rotate_left(20));
    let (p1, p2) = match torsion_basis(&curve, order_m, ell, 48, &mut rng) {
        Torsion::Rank2(p1, p2) => (p1, p2),
        Torsion::Rank1 => return Err("vertex is on the ℓ-floor: no descending edge".into()),
        Torsion::None => return Err("ℓ does not divide #E(F_{q^m})".into()),
        Torsion::Undecided => return Err("ℓ-torsion rank undecided within budget".into()),
    };
    let gens = kernel_generators(&curve, &p1, &p2, ell);
    let mut loops = 0usize;
    for (tried, xg) in gens.iter().enumerate() {
        let ker = kernel_of(&curve, xg, ell);
        if !ker.t_is_rational {
            return Err("kernel trace left F_q: not a rational kernel".into());
        }
        let tv = ker.t[0];
        let a6p = a6 ^ tv ^ inst.gf.mul(tv, tv);
        if avoid.contains(&a6p) {
            loops += 1;
            continue;
        }
        let (xi, rational) = velu_x_from_abscissae(ext, &ker.abscissae, x_carry)
            .ok_or("carried abscissa is a kernel abscissa")?;
        if !rational {
            return Err("image abscissa left F_q".into());
        }
        // Where the volcano has more than one level, the codomain's rank
        // says whether this edge went down (rank stays 2 above the floor,
        // 1 at the floor).  Only a floor landing is demanded when asked.
        let cod = ExtCurve::new(ext, inst.a, a6p);
        let mut rng2 = StdRng::seed_from_u64(seed ^ a6p ^ 0xF100);
        let at_floor = match torsion_basis(&cod, order_m, ell, 48, &mut rng2) {
            Torsion::Rank1 => true,
            Torsion::Rank2(..) => false,
            _ => return Err("codomain rank undecided".into()),
        };
        if let Some(want) = want_floor {
            if want != at_floor {
                loops += 1;
                continue;
            }
        }
        return Ok((
            BinaryEdge {
                a6_from: a6,
                a6_to: a6p,
                ell,
                m: ext.m,
                x_image: xi,
                kernels_tried: tried + 1,
                loops_skipped: loops,
            },
            at_floor,
        ));
    }
    Err("every kernel was a loop or a backtrack".into())
}

fn edge_json(e: &BinaryEdge) -> Value {
    json!({"ell": e.ell, "ext_degree": e.m, "a6_from": e.a6_from, "a6_to": e.a6_to,
           "x_image": e.x_image, "kernels_tried": e.kernels_tried, "loops_skipped": e.loops_skipped})
}

/// Verify `#E′ = #E` and lift the carried abscissa to a point of order `r`.
fn certify_binary_vertex(
    inst: &BinaryInstance,
    a6: u64,
    x: u64,
    count_limit: u32,
    rng: &mut StdRng,
) -> Result<(u64, u64), String> {
    let n = inst.n;
    let order = inst.group_order;
    if n <= count_limit {
        let counted = count_order(&inst.irreducible, inst.a, a6);
        if counted != order {
            return Err(format!("#E′ = {counted} ≠ #E = {order}"));
        }
    }
    let curve = BinCurve::with_order(n, &inst.irreducible, inst.a as u8, a6, order)
        .ok_or("codomain rejected by the curve constructor")?;
    let g = curve.group();
    let mut ops = GroupOps::default();
    for _ in 0..8 {
        let xr = rng.gen_range(1..(1u64 << n));
        if let Some(&pt) = curve.points_with_x(xr).first() {
            if !g.mul(&mut ops, pt, order).infinity {
                return Err("[#E] P ≠ O on the codomain".into());
            }
        }
    }
    lift_order_r(&curve, x, inst.r).ok_or("image abscissa carries no point of order r".into())
}

fn koblitz_ladder(a: u8, n: u32, seed: u64) -> Result<Value, String> {
    let started = Instant::now();
    let inst = koblitz_instance(a, n).ok_or("no Koblitz instance")?;
    let q = 1u64 << n;
    let order = inst.group_order;
    let t = (q as i64 + 1) - order as i64;
    let disc = t * t - 4 * q as i64;
    let (fund, cond) = fundamental_discriminant_and_conductor(disc);
    if fund != -7 {
        return Err(format!("Koblitz discriminant is not −7·f²: {fund}"));
    }
    let factors = factor_u64(cond as u64);
    let mut rng = StdRng::seed_from_u64(seed ^ n as u64);
    let mut a6 = inst.b;
    let mut x = inst.generator.x;
    let mut edges = Vec::new();
    let mut gap = 1u64;
    let mut levels = vec![json!({"role": "crater", "conductor": 1u64, "a6": a6,
        "end_evidence": "defined over F_2: End = Z[τ] = O_K",
        "spec": json!({"kind": "koblitz", "a": a, "n": n}),
        "explicit": binary_explicit_json(n, &inst.irreducible, inst.a, a6, order, inst.r, (inst.generator.x, inst.generator.y))})];
    let mut floor_certs = Vec::new();
    for &(ell, e) in &factors {
        if e == 0 {
            return Err(format!(
                "conductor factor {ell} is composite beyond the trial bound"
            ));
        }
        let lam =
            scalar_eigenvalue(q, t, ell).ok_or(format!("no scalar eigenvalue at ℓ = {ell}"))?;
        let m = mult_order(lam, ell) as usize;
        let ext = Ext::new(&inst.irreducible, m, seed ^ ell);
        let order_m = order_over_extension(q, t, m);
        let mut prev: Vec<u64> = vec![a6];
        for step in 0..e {
            let last = step + 1 == e;
            let (edge, at_floor) =
                descend_binary(&inst, &ext, ell, &order_m, a6, x, &prev, Some(last), seed)?;
            let g = certify_binary_vertex(&inst, edge.a6_to, edge.x_image, 25, &mut rng)?;
            prev.push(edge.a6_to);
            a6 = edge.a6_to;
            x = edge.x_image;
            gap *= ell;
            edges.push(edge_json(&edge));
            if last {
                floor_certs.push(json!({"ell": ell, "ext_degree": m, "rank": 1}));
            }
            levels.push(
                json!({"role": "level", "conductor": gap, "a6": a6, "ell": ell, "step": step + 1,
                "at_ell_floor": at_floor,
                "spec": binary_explicit_json(n, &inst.irreducible, inst.a, a6, order, inst.r, g)}),
            );
        }
    }
    let last = levels.len() - 1;
    if last == 0 {
        return Err("conductor 1: the class is a single vertex".into());
    }
    levels[last]["role"] = json!("floor");
    levels[last]["end_evidence"] = json!("rank-1 E[ℓ] over F_{q^m} for every ℓ | f_π: End = Z[π]");
    Ok(json!({
        "family": "koblitz", "a": a, "n": n, "trace": t, "group_order": order, "r": inst.r, "cofactor": inst.cofactor,
        "frobenius_conductor": cond, "conductor_factors": factors, "gap": gap,
        "edges": edges, "floor_certificates": floor_certs, "curves": levels,
        "seconds": started.elapsed().as_secs_f64(),
    }))
}

fn binary_random_pair(
    n: u32,
    seeds: u64,
    min_ell: u64,
    max_ell: u64,
    max_ext_bits: usize,
    seed: u64,
) -> Result<Value, String> {
    let started = Instant::now();
    let q = 1u64 << n;
    // Scan: the class whose conductor has the largest admissible prime ℓ.
    let mut best: Option<(u64, u64, usize, BinaryInstance)> = None;
    let mut scanned = 0u64;
    let mut koblitz_class_skipped = 0u64;
    // Traces of K_0 and K_1 over F_{2^n}: a random curve with one of them
    // (or its negative, the twist) lies in a Koblitz class, which the
    // Koblitz family already covers.
    let koblitz_traces: Vec<i64> = [-1i64, 1]
        .iter()
        .map(|&t1| {
            let (mut s0, mut s1) = (2i64, t1);
            for _ in 1..n {
                let s2 = t1 * s1 - 2 * s0;
                s0 = s1;
                s1 = s2;
            }
            s1.abs()
        })
        .collect();
    for s in 1..=seeds {
        let Some(inst) = random_binary_instance(n, s, 16) else {
            continue;
        };
        scanned += 1;
        let t = (q as i64 + 1) - inst.group_order as i64;
        if koblitz_traces.contains(&t.abs()) {
            koblitz_class_skipped += 1;
            continue;
        }
        let (_, cond) = fundamental_discriminant_and_conductor(t * t - 4 * q as i64);
        for (ell, e) in factor_u64(cond as u64) {
            if e != 1 || !(min_ell..=max_ell).contains(&ell) {
                continue;
            }
            let Some(lam) = scalar_eigenvalue(q, t, ell) else {
                continue;
            };
            let m = mult_order(lam, ell) as usize;
            if m * n as usize > max_ext_bits {
                continue;
            }
            if best.as_ref().is_none_or(|b| ell > b.1) {
                best = Some((s, ell, m, inst.clone_shallow()));
            }
        }
    }
    let (s, ell, m, inst) = best.ok_or(format!(
        "no class with an admissible ℓ among {scanned} curves"
    ))?;
    let t = (q as i64 + 1) - inst.group_order as i64;
    let (fund, cond) = fundamental_discriminant_and_conductor(t * t - 4 * q as i64);
    let ext = Ext::new(&inst.irreducible, m, seed ^ ell);
    let order_m = order_over_extension(q, t, m);
    let curve = ExtCurve::new(&ext, inst.a, inst.b);
    let mut rng = StdRng::seed_from_u64(seed ^ s);
    let base_spec = binary_explicit_json(
        n,
        &inst.irreducible,
        inst.a,
        inst.b,
        inst.group_order,
        inst.r,
        (inst.generator.x, inst.generator.y),
    );
    let (crater, floor, edge, direction) = match torsion_basis(&curve, &order_m, ell, 48, &mut rng)
    {
        Torsion::Rank2(..) => {
            let (edge, at_floor) = descend_binary(
                &inst,
                &ext,
                ell,
                &order_m,
                inst.b,
                inst.generator.x,
                &[inst.b],
                Some(true),
                seed,
            )?;
            if !at_floor {
                return Err("descent did not reach the floor".into());
            }
            let g = certify_binary_vertex(&inst, edge.a6_to, edge.x_image, 25, &mut rng)?;
            let floor_spec = binary_explicit_json(
                n,
                &inst.irreducible,
                inst.a,
                edge.a6_to,
                inst.group_order,
                inst.r,
                g,
            );
            (
                base_spec.clone(),
                floor_spec,
                edge,
                "found vertex is the crater; descended",
            )
        }
        Torsion::Rank1 => {
            // The unique rational line is the ascending kernel.
            let l = BigUint::from(ell);
            let mut v = 0u32;
            let mut u = order_m.clone();
            while (&u % &l).is_zero() {
                u /= &l;
                v += 1;
            }
            let _ = v;
            let qpt = loop {
                let xr = curve.random_x(&mut rng);
                let Some(yr) = curve.lift_x(&xr, &mut rng) else {
                    continue;
                };
                let mut qq = curve.mul_pt(&Pt::A(xr, yr), &u);
                if qq == Pt::O {
                    continue;
                }
                loop {
                    let next = curve.mul_pt(&qq, &l);
                    if next == Pt::O {
                        break;
                    }
                    qq = next;
                }
                break qq;
            };
            let Pt::A(xg, _) = qpt else { unreachable!() };
            let ker = kernel_of(&curve, &xg, ell);
            if !ker.t_is_rational {
                return Err("ascending kernel trace left F_q".into());
            }
            let tv = ker.t[0];
            let a6p = inst.b ^ tv ^ inst.gf.mul(tv, tv);
            let (xi, rational) = velu_x_from_abscissae(&ext, &ker.abscissae, inst.generator.x)
                .ok_or("generator in the kernel")?;
            if !rational {
                return Err("image abscissa left F_q".into());
            }
            let cod = ExtCurve::new(&ext, inst.a, a6p);
            let mut rng2 = StdRng::seed_from_u64(seed ^ a6p);
            if !matches!(
                torsion_basis(&cod, &order_m, ell, 48, &mut rng2),
                Torsion::Rank2(..)
            ) {
                return Err("ascended vertex is not of rank 2".into());
            }
            let g = certify_binary_vertex(&inst, a6p, xi, 25, &mut rng)?;
            let crater_spec = binary_explicit_json(
                n,
                &inst.irreducible,
                inst.a,
                a6p,
                inst.group_order,
                inst.r,
                g,
            );
            let edge = BinaryEdge {
                a6_from: inst.b,
                a6_to: a6p,
                ell,
                m,
                x_image: xi,
                kernels_tried: 1,
                loops_skipped: 0,
            };
            (
                crater_spec,
                base_spec.clone(),
                edge,
                "found vertex is on the floor; ascended",
            )
        }
        Torsion::None => return Err("ℓ ∤ #E(F_{q^m})".into()),
        Torsion::Undecided => return Err("rank undecided".into()),
    };
    Ok(json!({
        "family": "binary", "n": n, "seed_found": s, "a": inst.a, "b_found": inst.b, "trace": t,
        "group_order": inst.group_order, "r": inst.r, "cofactor": inst.cofactor,
        "fundamental_discriminant": fund, "frobenius_conductor": cond, "ell": ell, "ext_degree": m, "gap": ell,
        "direction": direction, "edge": edge_json(&edge),
        "curves": [
            {"role": "crater", "conductor": cond as u64 / ell, "end_evidence": "rank-2 E[ℓ] over F_{q^m}: ℓ ∤ [End : Z[π]] fails, so ℓ ∤ f_E", "spec": crater},
            {"role": "floor", "conductor": cond, "end_evidence": "rank-1 E[ℓ] over F_{q^m}: ℓ | f_E", "spec": floor},
        ],
        "scanned": scanned, "koblitz_class_skipped": koblitz_class_skipped, "seconds": started.elapsed().as_secs_f64(),
    }))
}

trait CloneShallow {
    fn clone_shallow(&self) -> Self;
}
impl CloneShallow for BinaryInstance {
    fn clone_shallow(&self) -> Self {
        random_binary_instance_again(self)
    }
}
fn random_binary_instance_again(inst: &BinaryInstance) -> BinaryInstance {
    // BinaryInstance is not Clone; rebuild it from its explicit form.
    crypto_lib::cryptanalysis::ecbench_large_prime::binary_explicit_instance(
        inst.n,
        modulus_mask(&inst.irreducible),
        inst.a,
        inst.b,
        inst.group_order,
        inst.r,
        inst.generator.x,
        inst.generator.y,
    )
    .expect("an instance rebuilds from its own fields")
}

// ── ecbench specs from the pairs ────────────────────────────────────

fn load_pairs(path: &str) -> Vec<Value> {
    match std::fs::read_to_string(path) {
        Ok(text) => serde_json::from_str::<Value>(&text)
            .ok()
            .and_then(|v| v["pairs"].as_array().cloned())
            .unwrap_or_default(),
        Err(_) => Vec::new(),
    }
}

fn spec_json(
    label: &str,
    description: &str,
    curves: Vec<Value>,
    arms: Vec<Value>,
    targets: u64,
    rounds: u64,
    timeout: u64,
) -> Value {
    json!({
        "schema": "ecbench.spec/v1",
        "label": label,
        "description": description,
        "workloads": {"curves": curves, "targets_per_curve": targets, "target_seed": 20261010},
        "arms": arms,
        "measurement": {"rounds": rounds, "warmup": 1, "order": "alternate", "seed": 1,
                        "isolation_required": "L2", "timeout_seconds": timeout}
    })
}

fn arm(name: &str, role: &str, id: &str, params: Value) -> Value {
    if params.is_null() {
        json!({"name": name, "role": role, "method": {"id": id}})
    } else {
        json!({"name": name, "role": role, "method": {"id": id, "params": params}})
    }
}

/// Binary factor-base dimension for a field degree: `2^{n − dim}`
/// algebraic solves per run (columns `≈ 2^{dim−1}`, hit rate
/// `≈ 2^{2·dim − n − 1}`), so 8, 9 and 11 at `n = 17, 19, 23` keep a run
/// to a few thousand solves.  Amendment to the protocol's "dimension 8",
/// recorded in the README.
fn binary_dimension(n: u64) -> u64 {
    match n {
        17 => 8,
        19 => 9,
        23 => 11,
        _ => n / 2 - 1,
    }
}

fn emit_specs(pairs_dir: &str, out: &str, targets: u64, rounds: u64) {
    std::fs::create_dir_all(out).expect("spec dir");
    // prime: crater, floor, isomorphic floor model, every pair
    let prime = load_pairs(&format!("{pairs_dir}/prime_pairs.json"));
    let mut by_bits: std::collections::BTreeMap<u64, Vec<Value>> =
        std::collections::BTreeMap::new();
    for pr in &prime {
        let entry = by_bits.entry(pr["bits"].as_u64().unwrap()).or_default();
        for c in pr["curves"].as_array().unwrap() {
            if c["role"] == "crater" || c["role"] == "floor" {
                entry.push(c["spec"].clone());
            }
        }
        entry.push(pr["floor_isomorphic_model"]["spec"].clone());
    }
    for (bits, curves) in by_bits {
        // The CM pairs have cofactors from 4 to 5408, and a fiber-blind
        // pair table holds only about 512/h sums inside the subgroup, so
        // the blind MITM exhausts on the large-cofactor pairs (the
        // interrupted prime_gap_22 session is the record).  The IC arm
        // across the gap is therefore the fiber-closed MITM, the same on
        // both curves of every pair.
        let prime_arms = vec![
            arm("rho-neg", "reference", "rho.negation", Value::Null),
            arm(
                "ic-mitm-closed",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": "prime-abscissa:size=32", "oracle": "mitm-fiber:fiber=closed,negation_folded=1"}),
            ),
            arm("rho-neg-aa", "control", "rho.negation", Value::Null),
        ];
        let spec = spec_json(
            &format!("isogenous pairs, prime fields at {bits} bits: crater, floor and an isomorphic floor model"),
            "Part B of research/isogeny_conductor_gap_ic_20261010/PROTOCOL.md: CM crater (End = O_K) and floor (End = Z[pi]) of every certified pair, plus an isomorphic model of the floor as the I-4 control; the same prime-abscissa fiber-closed MITM index calculus and the matched rho on every curve.",
            curves, prime_arms, targets, rounds, 900,
        );
        std::fs::write(
            format!("{out}/prime_gap_{bits}_fiber.json"),
            serde_json::to_string_pretty(&spec).unwrap(),
        )
        .unwrap();
    }

    // koblitz and binary: arms both ends can run
    let mut kob: Vec<Value> = load_pairs(&format!("{pairs_dir}/koblitz_k0_pairs.json"));
    kob.extend(load_pairs(&format!("{pairs_dir}/koblitz_k1_pairs.json")));
    let bin = load_pairs(&format!("{pairs_dir}/binary_pairs.json"));
    for (family, pairs) in [("koblitz", &kob), ("binary", &bin)] {
        // group by n so one spec carries one factor-base dimension
        let mut by_n: std::collections::BTreeMap<u64, Vec<Value>> =
            std::collections::BTreeMap::new();
        for pr in pairs {
            let n = pr["n"].as_u64().unwrap();
            let entry = by_n.entry(n).or_default();
            for c in pr["curves"].as_array().unwrap() {
                if c["role"] == "crater" || c["role"] == "floor" {
                    entry.push(c["spec"].clone());
                }
            }
        }
        for (n, curves) in by_n {
            let dim = binary_dimension(n);
            let arms = vec![
                arm("rho-neg", "reference", "rho.negation", Value::Null),
                arm(
                    "ic-f4",
                    "candidate",
                    "ic.pipeline",
                    json!({"factor_base": format!("binary-subspace:dimension={dim}"), "oracle": "descent-algebraic:m=2", "solver": "f4-f2"}),
                ),
                arm(
                    "ic-mitm",
                    "candidate",
                    "ic.pipeline",
                    json!({"factor_base": format!("binary-subspace:dimension={dim}"), "oracle": "mitm:negation_folded=1"}),
                ),
                arm("rho-neg-aa", "control", "rho.negation", Value::Null),
            ];
            let spec = spec_json(
                &format!("isogenous pairs, {family} n={n}: crater and floor under the same binary index calculus"),
                "Part B of research/isogeny_conductor_gap_ic_20261010/PROTOCOL.md: both ends of every certified pair under the subspace factor base with the F4 descent oracle (solving degree) and the MITM oracle (yield), and the generic rho on both.",
                curves, arms, targets, rounds, 1800,
            );
            std::fs::write(
                format!("{out}/{family}_gap_n{n}.json"),
                serde_json::to_string_pretty(&spec).unwrap(),
            )
            .unwrap();
        }
    }
    // koblitz craters only: the folds the floor cannot run
    let craters: Vec<Value> = kob
        .iter()
        .flat_map(|pr| {
            pr["curves"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|c| c["role"] == "crater")
                .map(|c| c["spec"].clone())
                .collect::<Vec<_>>()
        })
        .collect();
    if !craters.is_empty() {
        let arms = vec![
            arm(
                "rho-strong",
                "reference",
                "rho.signed_frobenius_strong",
                Value::Null,
            ),
            arm("rho-frob", "baseline", "rho.signed_frobenius", Value::Null),
            arm(
                "ic-orbit",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": "koblitz-orbit:divisor=0;1", "oracle": "mitm-frobenius:m=2"}),
            ),
        ];
        let spec = spec_json(
            "isogenous pairs, Koblitz craters: the Frobenius folds the floor cannot run",
            "Part B: on the crater K_a the signed-Frobenius rho and the tau-orbit factor base apply; the isogenous floor curve is not defined over F_2 and has neither, which is the one structural difference across the gap.",
            craters, arms, targets, rounds, 900,
        );
        std::fs::write(
            format!("{out}/koblitz_gap_crater_folds.json"),
            serde_json::to_string_pretty(&spec).unwrap(),
        )
        .unwrap();
    }
    eprintln!("specs written to {out}");
}

/// Part A: the fiber-aware relation-generation arms on curves with a
/// cofactor fiber, plus the prime-order control.
fn emit_fiber_specs(pairs_dir: &str, out: &str, targets: u64, rounds: u64) {
    std::fs::create_dir_all(out).expect("spec dir");
    // Binary (Koblitz and random): MITM arms, Frobenius arms (Koblitz
    // only) and algebraic arms, one spec per field degree.
    let koblitz: Vec<(u8, u32)> = vec![(1, 17), (0, 19), (1, 19), (0, 23)];
    for (a, n) in koblitz {
        let dim = binary_dimension(n as u64);
        let fb = format!("binary-subspace:dimension={dim}");
        let curves = vec![json!({"kind": "koblitz", "a": a, "n": n})];
        let arms = vec![
            arm("rho-neg", "reference", "rho.negation", Value::Null),
            arm(
                "mitm-blind",
                "baseline",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "mitm:negation_folded=1"}),
            ),
            arm(
                "mitm-lifts",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "mitm-fiber:fiber=lifts,negation_folded=1"}),
            ),
            arm(
                "mitm-closed",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "mitm-fiber:fiber=closed,negation_folded=1"}),
            ),
            arm(
                "frob-blind",
                "baseline",
                "ic.pipeline",
                json!({"factor_base": "koblitz-orbit:divisor=0;1", "oracle": "mitm-frobenius:m=2"}),
            ),
            arm(
                "frob-lifts",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": "koblitz-orbit:divisor=0;1", "oracle": "mitm-frobenius-fiber:fiber=lifts"}),
            ),
            arm(
                "frob-closed",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": "koblitz-orbit:divisor=0;1", "oracle": "mitm-frobenius-fiber:fiber=closed"}),
            ),
            arm(
                "f4-blind",
                "baseline",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "descent-algebraic:m=2", "solver": "f4-f2"}),
            ),
            arm(
                "f4-lifts",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "descent-algebraic-fiber:fiber=lifts", "solver": "f4-f2"}),
            ),
            arm("rho-neg-aa", "control", "rho.negation", Value::Null),
        ];
        let spec = spec_json(
            &format!("fiber-aware relation generation on K_{a} / F_2^{n}"),
            "Part A of research/isogeny_conductor_gap_ic_20261010/PROTOCOL.md: fiber-blind, all-lifts and fiber-closed MITM (plain and Frobenius-folded) and fiber-blind, all-lifts and fiber-combined F4 descent on one Koblitz curve; h = 4 (K_0) or 2 (K_1).",
            curves, arms, targets, rounds, 1800,
        );
        std::fs::write(
            format!("{out}/fiber_koblitz_k{a}_n{n}.json"),
            serde_json::to_string_pretty(&spec).unwrap(),
        )
        .unwrap();
    }
    for n in [17u32, 19] {
        let dim = binary_dimension(n as u64);
        let fb = format!("binary-subspace:dimension={dim}");
        let curves: Vec<Value> = (1..=4u64)
            .map(|seed| json!({"kind": "binary_random", "n": n, "seed": seed, "max_cofactor": 8}))
            .collect();
        let arms = vec![
            arm("rho-neg", "reference", "rho.negation", Value::Null),
            arm(
                "mitm-blind",
                "baseline",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "mitm:negation_folded=1"}),
            ),
            arm(
                "mitm-lifts",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "mitm-fiber:fiber=lifts,negation_folded=1"}),
            ),
            arm(
                "mitm-closed",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "mitm-fiber:fiber=closed,negation_folded=1"}),
            ),
            arm(
                "f4-blind",
                "baseline",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "descent-algebraic:m=2", "solver": "f4-f2"}),
            ),
            arm(
                "f4-lifts",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "descent-algebraic-fiber:fiber=lifts", "solver": "f4-f2"}),
            ),
            arm("rho-neg-aa", "control", "rho.negation", Value::Null),
        ];
        let spec = spec_json(
            &format!("fiber-aware relation generation on random binary curves, n = {n}"),
            "Part A: the same arms on four random binary curves with cofactor up to 8; the fiber size each oracle found is in its set-up counters.",
            curves, arms, targets, rounds, 1800,
        );
        std::fs::write(
            format!("{out}/fiber_binary_n{n}.json"),
            serde_json::to_string_pretty(&spec).unwrap(),
        )
        .unwrap();
    }
    // The combined algebraic system at n = 17 only: its F4 solve is two
    // orders of magnitude slower than the plain descent's (smoke run), so
    // it gets its own small spec rather than a place in every family's.
    {
        let dim = binary_dimension(17);
        let fb = format!("binary-subspace:dimension={dim}");
        let curves = vec![
            json!({"kind": "koblitz", "a": 1, "n": 17}),
            json!({"kind": "binary_random", "n": 17, "seed": 1, "max_cofactor": 8}),
            json!({"kind": "binary_random", "n": 17, "seed": 2, "max_cofactor": 8}),
        ];
        let arms = vec![
            arm("rho-neg", "reference", "rho.negation", Value::Null),
            arm(
                "f4-blind",
                "baseline",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "descent-algebraic:m=2", "solver": "f4-f2"}),
            ),
            arm(
                "f4-lifts",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "descent-algebraic-fiber:fiber=lifts", "solver": "f4-f2"}),
            ),
            arm(
                "f4-combined",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": fb, "oracle": "descent-algebraic-fiber:fiber=combined", "solver": "f4-f2"}),
            ),
        ];
        let spec = spec_json(
            "fiber-combined algebraic descent at n = 17",
            "Part A, arm A3: one Weil-descended system per target with the target abscissa free and the fiber polynomial adjoined, against the blind and all-lifts descents, on K_1 and two random binary curves over F_2^17.",
            curves, arms, 2.min(targets), 2.min(rounds), 3600,
        );
        std::fs::write(
            format!("{out}/fiber_combined_n17.json"),
            serde_json::to_string_pretty(&spec).unwrap(),
        )
        .unwrap();
    }
    // Prime: the h = 1 control and CM craters/floors with h = ℓ²·c.
    let prime = load_pairs(&format!("{pairs_dir}/prime_pairs.json"));
    let mut curves = vec![json!({"kind": "prime_search", "bits": 20, "seed": 59297})];
    for pr in &prime {
        if pr["bits"].as_u64() != Some(22) {
            continue;
        }
        let cof = pr["cofactor"].as_u64().unwrap_or(0);
        if !(2..=50).contains(&cof) {
            continue;
        }
        for c in pr["curves"].as_array().unwrap() {
            if c["role"] == "crater" || c["role"] == "floor" {
                curves.push(c["spec"].clone());
            }
        }
    }
    // The large-cofactor CM curves (h > 50): the blind and all-lifts arms
    // are capped at two million trials, so an arm that cannot find enough
    // in-subgroup pair sums exhausts within the session instead of
    // running for a minute per target.
    let mut large = Vec::new();
    for pr in &prime {
        if pr["bits"].as_u64() != Some(22) || pr["cofactor"].as_u64().unwrap_or(0) <= 50 {
            continue;
        }
        for c in pr["curves"].as_array().unwrap() {
            if c["role"] == "crater" || c["role"] == "floor" {
                large.push(c["spec"].clone());
            }
        }
    }
    if !large.is_empty() {
        let arms = vec![
            arm("rho-neg", "reference", "rho.negation", Value::Null),
            arm(
                "mitm-blind",
                "baseline",
                "ic.pipeline",
                json!({"factor_base": "prime-abscissa:size=32", "oracle": "mitm:negation_folded=1", "max_trials": "2000000"}),
            ),
            arm(
                "mitm-lifts",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": "prime-abscissa:size=32", "oracle": "mitm-fiber:fiber=lifts,negation_folded=1", "max_trials": "2000000"}),
            ),
            arm(
                "mitm-closed",
                "candidate",
                "ic.pipeline",
                json!({"factor_base": "prime-abscissa:size=32", "oracle": "mitm-fiber:fiber=closed,negation_folded=1"}),
            ),
            arm("rho-neg-aa", "control", "rho.negation", Value::Null),
        ];
        let spec = spec_json(
            "fiber-aware relation generation on prime CM curves with cofactor above 50",
            "Part A on the certified craters and floors at 22 bits whose cofactor exceeds 50: the fiber-blind and all-lifts MITM capped at two million trials against the fiber-closed MITM, which needs no cap.",
            large, arms, targets, rounds, 900,
        );
        std::fs::write(
            format!("{out}/fiber_prime_large_h.json"),
            serde_json::to_string_pretty(&spec).unwrap(),
        )
        .unwrap();
    }
    let arms = vec![
        arm("rho-neg", "reference", "rho.negation", Value::Null),
        arm(
            "mitm-blind",
            "baseline",
            "ic.pipeline",
            json!({"factor_base": "prime-abscissa:size=32", "oracle": "mitm:negation_folded=1"}),
        ),
        arm(
            "mitm-lifts",
            "candidate",
            "ic.pipeline",
            json!({"factor_base": "prime-abscissa:size=32", "oracle": "mitm-fiber:fiber=lifts,negation_folded=1"}),
        ),
        arm(
            "mitm-closed",
            "candidate",
            "ic.pipeline",
            json!({"factor_base": "prime-abscissa:size=32", "oracle": "mitm-fiber:fiber=closed,negation_folded=1"}),
        ),
        arm("rho-neg-aa", "control", "rho.negation", Value::Null),
    ];
    let spec = spec_json(
        "fiber-aware relation generation on prime curves: the h = 1 control and CM curves with h = ell^2 c",
        "Part A: fiber-blind, all-lifts and fiber-closed MITM on a prime-order curve (every mode must coincide) and on the certified CM craters and floors at 22 bits whose cofactor is between 2 and 50.",
        curves, arms, targets, rounds, 900,
    );
    std::fs::write(
        format!("{out}/fiber_prime.json"),
        serde_json::to_string_pretty(&spec).unwrap(),
    )
    .unwrap();
    eprintln!("fiber specs written to {out}");
}

/// The registry source for every explicit curve of the round
/// (`docs/curves/sources/*.json`, the `ecbench.curve-source/v1` layout
/// `scripts/build_curve_registry.py` reads).
fn emit_registry_source(pairs_dir: &str, out_file: &str) {
    let mut curves = Vec::new();
    let mut push = |name: String, spec: &Value| match spec["kind"].as_str().unwrap_or("") {
        "prime_explicit" => curves.push(json!({
            "name": name, "p": spec["p"], "a": spec["a"], "b": spec["b"],
            "group_order": spec["group_order"], "subgroup_order": spec["r"],
            "cofactor": spec["group_order"].as_u64().unwrap() / spec["r"].as_u64().unwrap(),
            "generator": {"x": spec["gx"], "y": spec["gy"]},
        })),
        "binary_explicit" => {
            let n = spec["n"].as_u64().unwrap();
            let modulus = spec["modulus"].as_u64().unwrap();
            let low: Vec<u64> = (0..n).filter(|k| (modulus >> k) & 1 == 1).collect();
            curves.push(json!({
                "name": name,
                "field": {"kind": "binary", "degree": n, "polynomial_low_terms": low},
                "a": spec["a"], "b": spec["b"],
                "group_order": spec["group_order"], "subgroup_order": spec["r"],
                "cofactor": spec["group_order"].as_u64().unwrap() / spec["r"].as_u64().unwrap(),
                "generator": {"x": spec["gx"], "y": spec["gy"]},
            }))
        }
        _ => {}
    };
    for file in [
        "prime_pairs.json",
        "koblitz_k0_pairs.json",
        "koblitz_k1_pairs.json",
        "binary_pairs.json",
    ] {
        let pairs = load_pairs(&format!("{pairs_dir}/{file}"));
        for pr in &pairs {
            let family = pr["family"].as_str().unwrap_or("?");
            let tag = match family {
                "prime" => format!(
                    "prime-D{}-l{}-h{}-{}b",
                    pr["d_k"].as_i64().unwrap().abs(),
                    pr["ell"],
                    pr["h"],
                    pr["bits"]
                ),
                "koblitz" => format!("koblitz-k{}-n{}", pr["a"], pr["n"]),
                _ => format!("binary-n{}-l{}", pr["n"], pr["ell"]),
            };
            for c in pr["curves"].as_array().unwrap() {
                let role = c["role"].as_str().unwrap_or("level");
                let lvl = c["level"]
                    .as_u64()
                    .map(|l| format!("-L{l}"))
                    .unwrap_or_default();
                push(format!("isogeny-gap-{tag}-{role}{lvl}"), &c["spec"]);
            }
            if pr["floor_isomorphic_model"].is_object() {
                push(
                    format!("isogeny-gap-{tag}-floor-iso"),
                    &pr["floor_isomorphic_model"]["spec"],
                );
            }
        }
    }
    let doc = json!({
        "schema": "ecbench.curve-source/v1",
        "source": "research/isogeny_conductor_gap_ic_20261010/pairs (examples/isogeny_gap_pairs.rs)",
        "curves": curves,
    });
    std::fs::write(out_file, serde_json::to_string_pretty(&doc).unwrap())
        .expect("write registry source");
    eprintln!("registry source written to {out_file}");
}

// ── main ────────────────────────────────────────────────────────────

fn arg(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1).cloned())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let family = arg(&args, "--family").unwrap_or_else(|| "prime".into());
    let out = arg(&args, "--out").unwrap_or_else(|| ".".into());
    let out_name = arg(&args, "--out-name");
    let seed: u64 = arg(&args, "--seed")
        .and_then(|v| v.parse().ok())
        .unwrap_or(20261010);
    std::fs::create_dir_all(&out).expect("output directory");
    let mut results = Vec::new();
    let mut failures = Vec::new();
    match family.as_str() {
        "prime" => {
            let bits: Vec<u32> = arg(&args, "--bits")
                .map(|v| v.split(',').filter_map(|s| s.parse().ok()).collect())
                .unwrap_or_else(|| vec![22, 26]);
            // (D_K, ℓ, h): the gap is ℓ^h.
            let configs: Vec<(i64, u64, u32)> = vec![
                (-7, 3, 6),
                (-7, 5, 4),
                (-7, 7, 3),
                (-7, 2, 8),
                (-7, 13, 2),
                (-8, 3, 6),
                (-8, 5, 4),
                (-8, 11, 2),
                (-11, 3, 6),
                (-11, 7, 3),
                (-11, 2, 8),
                (-19, 3, 6),
                (-19, 5, 4),
                (-19, 11, 2),
                (-3, 5, 4),
                (-3, 7, 3),
                (-3, 2, 8),
                (-4, 3, 6),
                (-4, 7, 3),
                (-4, 5, 4),
            ];
            for &b in &bits {
                for &(dk, ell, h) in &configs {
                    let label = format!("D={dk} ell={ell} h={h} bits={b}");
                    eprintln!("prime: {label}");
                    match prime_ladder(dk, ell, h, b, seed) {
                        Ok(v) => {
                            eprintln!(
                                "  ok: p={} gap={} in {:.1}s",
                                v["p"],
                                v["gap"],
                                v["seconds"].as_f64().unwrap()
                            );
                            results.push(v);
                        }
                        Err(e) => {
                            eprintln!("  FAILED: {e}");
                            failures.push(json!({"family": "prime", "d_k": dk, "ell": ell, "h": h, "bits": b, "error": e}));
                        }
                    }
                }
            }
        }
        "koblitz" => {
            let ns: Vec<u32> = arg(&args, "--n")
                .map(|v| v.split(',').filter_map(|s| s.parse().ok()).collect())
                .unwrap_or_else(|| vec![19, 23, 17, 20]);
            let a: u8 = arg(&args, "--a").and_then(|v| v.parse().ok()).unwrap_or(0);
            for n in ns {
                eprintln!("koblitz: K_{a} n={n}");
                match koblitz_ladder(a, n, seed) {
                    Ok(v) => {
                        eprintln!(
                            "  ok: gap={} in {:.1}s",
                            v["gap"],
                            v["seconds"].as_f64().unwrap()
                        );
                        results.push(v);
                    }
                    Err(e) => {
                        eprintln!("  FAILED: {e}");
                        failures.push(json!({"family": "koblitz", "a": a, "n": n, "error": e}));
                    }
                }
            }
        }
        "binary" => {
            let ns: Vec<u32> = arg(&args, "--n")
                .map(|v| v.split(',').filter_map(|s| s.parse().ok()).collect())
                .unwrap_or_else(|| vec![17, 19, 23]);
            let seeds: u64 = arg(&args, "--seeds")
                .and_then(|v| v.parse().ok())
                .unwrap_or(400);
            let min_ell: u64 = arg(&args, "--min-ell")
                .and_then(|v| v.parse().ok())
                .unwrap_or(11);
            let max_ell: u64 = arg(&args, "--max-ell")
                .and_then(|v| v.parse().ok())
                .unwrap_or(1200);
            let max_ext: usize = arg(&args, "--max-ext-bits")
                .and_then(|v| v.parse().ok())
                .unwrap_or(3000);
            for n in ns {
                eprintln!("binary: n={n}");
                match binary_random_pair(n, seeds, min_ell, max_ell, max_ext, seed) {
                    Ok(v) => {
                        eprintln!(
                            "  ok: ell={} m={} in {:.1}s",
                            v["ell"],
                            v["ext_degree"],
                            v["seconds"].as_f64().unwrap()
                        );
                        results.push(v);
                    }
                    Err(e) => {
                        eprintln!("  FAILED: {e}");
                        failures.push(json!({"family": "binary", "n": n, "error": e}));
                    }
                }
            }
        }
        "registry-source" => {
            let pairs_dir = arg(&args, "--pairs-dir").unwrap_or_else(|| out.clone());
            let file = arg(&args, "--out-file").expect("--out-file FILE");
            emit_registry_source(&pairs_dir, &file);
            return;
        }
        "fiber-specs" => {
            let pairs_dir = arg(&args, "--pairs-dir").unwrap_or_else(|| out.clone());
            let targets: u64 = arg(&args, "--targets")
                .and_then(|v| v.parse().ok())
                .unwrap_or(4);
            let rounds: u64 = arg(&args, "--rounds")
                .and_then(|v| v.parse().ok())
                .unwrap_or(3);
            emit_fiber_specs(&pairs_dir, &out, targets, rounds);
            return;
        }
        "specs" => {
            let pairs_dir = arg(&args, "--pairs-dir").unwrap_or_else(|| out.clone());
            let targets: u64 = arg(&args, "--targets")
                .and_then(|v| v.parse().ok())
                .unwrap_or(4);
            let rounds: u64 = arg(&args, "--rounds")
                .and_then(|v| v.parse().ok())
                .unwrap_or(3);
            emit_specs(&pairs_dir, &out, targets, rounds);
            return;
        }
        other => {
            eprintln!("unknown family {other}");
            std::process::exit(2);
        }
    }
    let path = format!(
        "{out}/{}.json",
        out_name.unwrap_or_else(|| format!("{family}_pairs"))
    );
    std::fs::write(
        &path,
        serde_json::to_string_pretty(
            &json!({"schema": "isogeny_gap_pairs/v1", "family": family, "seed": seed,
            "pairs": results, "failures": failures}),
        )
        .unwrap(),
    )
    .expect("write");
    eprintln!("wrote {path}");
}
