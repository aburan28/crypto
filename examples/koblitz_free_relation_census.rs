//! Free-relation census of the compact-orbit pair-root table.
//!
//! Lead B of `RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md` §3: two
//! pair-states with the same x-coordinate up to Frobenius are a 4-term
//! relation `P_i ± π^k P_j = ±π^t (P_u ± π^v P_w)` among the factor-base
//! orbit representatives, found by one sort and no probes. The note predicts
//! `≈ n⁴K⁴/r` of them; counting Frobenius classes of subgroup x-values
//! (`≈ r/(2n)` bins) instead gives `≈ K⁴n³/r`. This example counts them.
//!
//! It is written from the group law alone (library `FastBinaryCurve` /
//! `FastBinaryCurve128` and their fields), not from the IC pipeline's
//! normal-basis root solver, so it is an independent census of the same
//! table `examples/koblitz_orbit_dlp_fast.rs` builds:
//!
//! 1. for every unordered state `(i ≤ j, k)` — `k ≤ (n − 1)/2` when `i = j`,
//!    which removes the `(i, j, k) ~ (j, i, −k)` duplicates — form
//!    `S± = P_i ± π^k P_j` and key it by the least Frobenius rotation of
//!    `x(S±)`;
//! 2. sort; every key shared by two states is a collision;
//! 3. for each collision find `t, η` with `S₁ = η·π^t S₂` (a full point
//!    comparison: this *is* the group identity), write the row
//!    `e_i + ε λ^k e_j − η λ^t (e_u + ε' λ^v e_w)` over `Z/r`, drop rows that
//!    vanish mod `r`, and check each kept row as `Σ [c]P = O` by scalar
//!    multiplication;
//! 4. rank of the kept rows mod `r`.
//!
//! Usage: `koblitz_free_relation_census <base_header.jsonl> <out.json>`, or
//! `construct:<n>:<a>:<K>:<seed>` for a synthetic test base of `K` seeded
//! subgroup points at `n ≤ 63` (self-check of the census logic only).
//! (a `point_defined_factor_base` header, as committed for the n = 61, 73 and
//! 83 single-target fixtures).
use crypto_lib::cryptanalysis::koblitz_fast_arith::{FastBinaryCurve, FastBinaryCurve128};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};
use rayon::prelude::*;
use serde_json::{json, Value};
use std::time::Instant;

type P = Option<(u128, u128)>;

/// The two word widths behind one interface; points are carried as `u128`.
enum Arith {
    Narrow(FastBinaryCurve, Gf2),
    Wide(FastBinaryCurve128),
}

impl Arith {
    fn add(&self, p: P, q: P) -> P {
        match self {
            Arith::Narrow(c, _) => c
                .add(
                    p.map(|(x, y)| (x as u64, y as u64)),
                    q.map(|(x, y)| (x as u64, y as u64)),
                )
                .map(|(x, y)| (x as u128, y as u128)),
            Arith::Wide(c) => c.add(p, q),
        }
    }

    fn sqr(&self, a: u128) -> u128 {
        match self {
            Arith::Narrow(_, gf) => gf.sqr(a as u64) as u128,
            Arith::Wide(c) => c.gf.sqr(a),
        }
    }

    fn frob(&self, p: P, k: usize) -> P {
        p.map(|(mut x, mut y)| {
            for _ in 0..k {
                x = self.sqr(x);
                y = self.sqr(y);
            }
            (x, y)
        })
    }

    fn neg(p: P) -> P {
        p.map(|(x, y)| (x, x ^ y))
    }

    fn mul(&self, p: P, k: &BigUint) -> P {
        match self {
            Arith::Narrow(c, _) => c
                .scalar_mul(p.map(|(x, y)| (x as u64, y as u64)), k)
                .map(|(x, y)| (x as u128, y as u128)),
            Arith::Wide(c) => c.scalar_mul(p, k),
        }
    }

    /// Least Frobenius rotation of `x`, and the rotation that reaches it.
    fn canonical(&self, x: u128, n: usize) -> (u128, u8) {
        let (mut best, mut at, mut v) = (x, 0u8, x);
        for s in 1..n {
            v = self.sqr(v);
            if v < best {
                best = v;
                at = s as u8;
            }
        }
        (best, at)
    }
}

fn mulmod(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 * b as u128) % r as u128) as u64
}

fn powmod(a: u64, mut e: u64, r: u64) -> u64 {
    let (mut acc, mut base) = (1u64, a % r);
    while e > 0 {
        if e & 1 == 1 {
            acc = mulmod(acc, base, r);
        }
        base = mulmod(base, base, r);
        e >>= 1;
    }
    acc
}

fn rank_mod(mut rows: Vec<Vec<u64>>, cols: usize, r: u64) -> usize {
    let mut rank = 0;
    for c in 0..cols {
        let Some(p) = (rank..rows.len()).find(|&i| rows[i][c] != 0) else {
            continue;
        };
        rows.swap(rank, p);
        let inv = powmod(rows[rank][c], r - 2, r);
        let pivot: Vec<u64> = rows[rank].iter().map(|&v| mulmod(v, inv, r)).collect();
        rows[rank] = pivot.clone();
        rows[rank + 1..].par_iter_mut().for_each(|row| {
            let f = row[c];
            if f != 0 {
                for (x, &y) in row.iter_mut().zip(&pivot) {
                    *x = (*x + r - mulmod(f, y, r)) % r;
                }
            }
        });
        rank += 1;
    }
    rank
}

fn word(v: &Value) -> u128 {
    match v {
        Value::String(s) => s.parse().expect("field value"),
        _ => v.as_u64().expect("field value") as u128,
    }
}

/// `(i, j, k, sign)` packed; `sign = 1` for `P_i − π^k P_j`.
fn pack(i: usize, j: usize, k: usize, minus: bool) -> u64 {
    i as u64 | (j as u64) << 16 | (k as u64) << 32 | (minus as u64) << 48
}

fn unpack(v: u64) -> (usize, usize, usize, bool) {
    (
        (v & 0xffff) as usize,
        (v >> 16 & 0xffff) as usize,
        (v >> 32 & 0xffff) as usize,
        v >> 48 & 1 == 1,
    )
}

/// A synthetic header: `K` seeded points of `⟨G⟩` (test mode, `n ≤ 63`).
fn construct(spec: &str) -> Value {
    let v: Vec<u64> = spec.split(':').map(|t| t.parse().unwrap()).collect();
    let (n, a, kk, mut seed) = (v[0] as u32, v[1] as u8, v[2] as usize, v[3]);
    let curve = KoblitzCurve::new(a, n).expect("admitted Koblitz rung");
    let fast = FastBinaryCurve::new(&curve.curve.irreducible, a as u64).unwrap();
    let gf = Gf2::new(&curve.curve.irreducible);
    let b = gf.from_element(&curve.curve.b);
    let mut reps = vec![];
    while reps.len() < kk {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let x = (seed >> 11) & ((1u64 << n) - 1);
        if let Some(&p) = fast.points_with_x(b, x).first() {
            if let Some((x, y)) = fast.scalar_mul(p, &curve.cofactor) {
                reps.push(json!([x, y]));
            }
        }
    }
    json!({"kind": "point_defined_factor_base", "n": n, "a": a,
           "subgroup_order": curve.subgroup_order.to_u64().unwrap(),
           "base_hash": format!("construct:{spec}"), "factor_base_representatives": reps})
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 3, "usage: <base_header.jsonl> <out.json>");
    let started = Instant::now();
    let header: Value = match args[1].strip_prefix("construct:") {
        Some(spec) => construct(spec),
        None => {
            let bytes = std::fs::read(&args[1]).expect("read base header");
            let line = bytes.split(|&b| b == b'\n').next().unwrap();
            serde_json::from_slice(line).expect("parse base header")
        }
    };
    assert_eq!(header["kind"], "point_defined_factor_base");
    let n = header["n"].as_u64().unwrap() as u32;
    let a = header["a"].as_u64().unwrap() as u8;
    let curve = KoblitzCurve::new(a, n).expect("admitted Koblitz rung");
    let r = curve.subgroup_order.to_u64().expect("r fits u64");
    assert_eq!(header["subgroup_order"].as_u64(), Some(r));
    let lambda = curve.lambda.to_u64().unwrap();
    let irr = &curve.curve.irreducible;
    let ar = if n <= 63 {
        Arith::Narrow(FastBinaryCurve::new(irr, a as u64).unwrap(), Gf2::new(irr))
    } else {
        Arith::Wide(FastBinaryCurve128::new(irr, a as u128).unwrap())
    };
    let reps: Vec<P> = header["factor_base_representatives"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| Some((word(&c[0]), word(&c[1]))))
        .collect();
    let kk = reps.len();
    let nn = n as usize;
    let rb = BigUint::from(r);
    let lb = BigUint::from(lambda);
    // Every representative is in ⟨G⟩ and π acts on it as λ.
    let in_subgroup = reps
        .par_iter()
        .all(|&p| ar.mul(p, &rb).is_none() && ar.frob(p, 1) == ar.mul(p, &lb));
    assert!(in_subgroup, "representatives must lie in the λ-eigenspace of order r");
    let shifted: Vec<Vec<P>> = reps
        .iter()
        .map(|&p| (0..nn).map(|k| ar.frob(p, k)).collect())
        .collect();

    // 1. states and keys
    let mut records: Vec<(u128, u64)> = (0..kk)
        .into_par_iter()
        .flat_map_iter(|i| {
            let mut out = Vec::new();
            for j in i..kk {
                let kmax = if i == j { (nn - 1) / 2 } else { nn - 1 };
                for k in 0..=kmax {
                    let q = shifted[j][k];
                    for minus in [false, true] {
                        if i == j && k == 0 && minus {
                            continue; // P − P = O
                        }
                        let s = ar.add(reps[i], if minus { Arith::neg(q) } else { q });
                        if let Some((x, _)) = s {
                            out.push((ar.canonical(x, nn).0, pack(i, j, k, minus)));
                        }
                    }
                }
            }
            out.into_iter()
        })
        .collect();
    let states = records.len();
    let t_states = started.elapsed().as_secs_f64();
    records.par_sort_unstable();

    // 2. collisions
    let mut groups: Vec<Vec<u64>> = vec![];
    let mut idx = 0;
    while idx < records.len() {
        let mut end = idx + 1;
        while end < records.len() && records[end].0 == records[idx].0 {
            end += 1;
        }
        if end - idx > 1 {
            groups.push(records[idx..end].iter().map(|&(_, m)| m).collect());
        }
        idx = end;
    }
    drop(records);
    let colliding_states: usize = groups.iter().map(|g| g.len()).sum();

    // 3. relations
    let point = |m: u64| -> P {
        let (i, j, k, minus) = unpack(m);
        let q = shifted[j][k];
        ar.add(reps[i], if minus { Arith::neg(q) } else { q })
    };
    let lam_pow: Vec<u64> = (0..nn as u64).map(|k| powmod(lambda, k, r)).collect();
    let term = |m: u64, scale: u64, row: &mut Vec<u64>| {
        let (i, j, k, minus) = unpack(m);
        row[i] = (row[i] + scale) % r;
        let c = mulmod(scale, lam_pow[k], r);
        row[j] = (row[j] + if minus { r - c } else { c }) % r;
    };
    let mut identity_failures = 0usize;
    let mut trivial = 0usize;
    let mut rows: Vec<Vec<u64>> = vec![];
    for g in &groups {
        let s1 = point(g[0]);
        for &m2 in &g[1..] {
            let s2 = point(m2);
            let found = (0..nn).find_map(|t| {
                let f = ar.frob(s2, t);
                if f == s1 {
                    Some((t, 1u64))
                } else if f == Arith::neg(s1) {
                    Some((t, r - 1))
                } else {
                    None
                }
            });
            let Some((t, eta)) = found else {
                identity_failures += 1;
                continue;
            };
            let mut row = vec![0u64; kk];
            term(g[0], 1, &mut row);
            term(m2, r - mulmod(eta, lam_pow[t], r), &mut row);
            if row.iter().all(|&v| v == 0) {
                trivial += 1;
            } else {
                rows.push(row);
            }
        }
    }
    // Each kept row as a group identity: Σ [c]P = O.
    let row_failures = rows
        .par_iter()
        .filter(|row| {
            let mut acc: P = None;
            for (c, &v) in row.iter().enumerate() {
                if v != 0 {
                    acc = ar.add(acc, ar.mul(reps[c], &BigUint::from(v)));
                }
            }
            acc.is_some()
        })
        .count();
    let nontrivial = rows.len();
    let rank = rank_mod(rows, kk, r);

    let (log_r, log_k, log_n) = ((r as f64).log2(), (kk as f64).log2(), (n as f64).log2());
    let out = json!({
        "schema": "koblitz_free_relation_census/v1",
        "base_header": args[1],
        "base_hash": header["base_hash"],
        "n": n, "a": a, "r": r, "log2_r": log_r, "lambda": lambda,
        "orbit_columns": kk,
        "states_with_affine_sum": states,
        "collision_groups": groups.len(),
        "states_in_collision_groups": colliding_states,
        "max_group_size": groups.iter().map(|g| g.len()).max().unwrap_or(0),
        "identity_failures": identity_failures,
        "trivial_rows": trivial,
        "nontrivial_relations": nontrivial,
        "row_check_failures": row_failures,
        "rank_mod_r": rank,
        "rank_needed": kk - 1,
        "prediction_note_n4K4_over_r": (4.0 * log_n + 4.0 * log_k - log_r).exp2(),
        "prediction_classes_K4n3_over_r": (3.0 * log_n + 4.0 * log_k - log_r).exp2(),
        "seconds_states": t_states,
        "seconds_total": started.elapsed().as_secs_f64(),
    });
    let text = serde_json::to_string_pretty(&out).unwrap();
    std::fs::write(&args[2], &text).expect("write out");
    println!("{text}");
    assert!(!r.is_zero());
}
