//! **First binary-regime sub-`2^{2ℓ}` decomposition oracle at `ℓ = 8`.**
//!
//! Ledger target (regime A, binary non-Koblitz, stage `decomposition`):
//! a sub-`2^{2ℓ}` oracle at fixed `ℓ = 8` on ≤ 64 targets with median
//! per-target cost at most half of pairs-and-solve, zero disagreements
//! against the group check, and the measurement-schema fields logged.
//!
//! Design.  The reference oracle is [`SubspaceOracle`] (Gaudry
//! pairs-and-solve): it enumerates all `2^ℓ(2^ℓ+1)/2` pairs of subspace
//! abscissae and solves a quartic for the third, `Θ(2^{2ℓ})` per target.
//! The candidate is a **group-law meet-in-the-middle pair table**: all
//! sums `Pᵢ + Pⱼ` of factor-base *points* are precomputed once
//! (target-independent, charged to the factor-base stage), and a query
//! walks the `|F| ≈ 2^ℓ` base points `P₃`, computing `R − P₃` and looking
//! it up.  Per-target work is `Θ(2^ℓ)` point additions — `2^ℓ` below the
//! `2^{2ℓ}` frontier — at `Θ(|F|²)` retained memory, exactly the trade the
//! ledger's required factor-base fields (`materialized`,
//! `construction_wall_ms`, `retained_bytes`) describe.
//!
//! Both arms decide the same question on the same targets: does the
//! target decompose over the factor base at the abscissa level (either
//! sign of the target point, which is what `S₄` decides)?  The MITM arm
//! answers it exactly — querying both `R` and `−R` — and every witness is
//! re-added in the group before it is reported, so `None` from both
//! signs is a proof (the table is complete).
//!
//! The curve is a genuinely non-Koblitz binary curve
//! `y² + xy = x³ + a·x² + b` with `a = z ∉ F₂`, `b = 1` over `F_{2^25}`
//! (`n` odd so the half-trace solves the quadratic for `y`); the base is
//! the low-order subspace `⟨1, z, …, z⁷⟩` of dimension `ℓ = 8`.
//!
//! ```bash
//! cargo run --release --example binary_mitm_oracle_l8
//! ```

use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::koblitz_fast::{FastCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use crypto_lib::cryptanalysis::semaev_decomp::{Gf2, SubspaceOracle};
use num_bigint::BigUint;
use num_traits::One;
use std::collections::HashMap;
use std::time::Instant;

/// Field degree.  Odd, so [`F2mElement::half_trace`] solves `t²+t=c`;
/// `3ℓ = 24 ≤ 25` keeps the decomposition probability well away from both
/// 0 and 1, so the 64 targets carry a genuine found/refuted mix.
const N: u32 = 25;
/// Factor-base subspace dimension (the fixed ledger rung).
const L: u32 = 8;
/// Curve constants: `a = z` (non-Koblitz, `a ∉ F₂`), `b = 1 ≠ 0`.
const A_BITS: u64 = 0b10;
const B: u64 = 1;
/// Number of targets: 48 natural + 16 planted-satisfiable.
const NATURAL_TARGETS: usize = 48;
const PLANTED_TARGETS: usize = 16;

// ── tiny deterministic RNG (same convention as semaev_decomp_bench) ──

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.0
    }
}

/// Solve `y² + x·y = c` for `y` on a binary curve (both roots, or none).
fn solve_y(x: &F2mElement, c: &F2mElement, irr: &IrreduciblePoly) -> Vec<F2mElement> {
    if x.is_zero() {
        // y² = c: y = c^{2^{n-1}} (square root in characteristic two).
        let y = c.square_k_times(N - 1, irr);
        return if y.square(irr) == *c && !y.is_zero() {
            vec![y]
        } else {
            vec![]
        };
    }
    let x2 = x.square(irr);
    let Some(_x_inv) = x.flt_inverse(irr) else {
        return vec![];
    };
    let rhs = c.mul(&x2.flt_inverse(irr).unwrap(), irr);
    if rhs.trace(irr) != 0 {
        return vec![];
    }
    let Some(t) = rhs.half_trace(irr) else {
        return vec![];
    };
    let y0 = x.mul(&t, irr);
    let y1 = y0.add(x);
    let mut out = vec![y0.clone()];
    if y1 != y0 {
        out.push(y1);
    }
    out.retain(|y| y.square(irr).add(&x.mul(y, irr)) == *c);
    out
}

/// Build the curve, its fast (one-word) twin, and the factor-base points
/// over the low-order `ℓ`-dimensional subspace.
fn setup() -> (FastCurve, Gf2, SubspaceOracle, Vec<FastPoint>, F2mElement) {
    let irr = find_irreducible_sparse(N).expect("sparse irreducible of degree 25");
    let a = F2mElement::from_biguint(&BigUint::from(A_BITS), N);
    let b_elt = F2mElement::from_biguint(&BigUint::from(B), N);

    // A valid on-curve point to satisfy the BinaryCurve struct; the fast
    // twin only reads m, the modulus, a and b.
    let mut generator = None;
    let mut probe = 1u64;
    while generator.is_none() {
        let x = F2mElement::from_biguint(&BigUint::from(probe), N);
        let c = curve_rhs(&x, &a, &b_elt, &irr);
        if let Some(y) = solve_y(&x, &c, &irr).first() {
            generator = Some(BinaryPoint::Affine {
                x: x.clone(),
                y: y.clone(),
            });
        }
        probe += 1;
    }
    let curve = BinaryCurve {
        m: N,
        irreducible: irr.clone(),
        a: a.clone(),
        b: b_elt.clone(),
        generator: generator.unwrap(),
        order: BigUint::one(),
        cofactor: BigUint::one(),
    };
    let fast = FastCurve::new(&curve).expect("n = 25 fits the packed path");
    let gf = Gf2::new(&irr);
    let oracle = SubspaceOracle::low_order(L, B, &gf);

    // Factor base: every curve point whose abscissa lies in the subspace.
    let mut base: Vec<FastPoint> = Vec::new();
    for v in 0..(1u64 << L) {
        let x = F2mElement::from_biguint(&BigUint::from(v), N);
        let c = curve_rhs(&x, &a, &b_elt, &irr);
        for y in solve_y(&x, &c, &irr) {
            let p = FastPoint::affine(gf.from_element(&x), gf.from_element(&y));
            assert!(fast.is_on_curve(p), "lifted base point must be on curve");
            base.push(p);
        }
    }

    (fast, gf, oracle, base, b_elt)
}

/// Right-hand side `x³ + a·x² + b` of the curve equation.
fn curve_rhs(x: &F2mElement, a: &F2mElement, b: &F2mElement, irr: &IrreduciblePoly) -> F2mElement {
    let x2 = x.square(irr);
    let x3 = x2.mul(x, irr);
    x3.add(&a.mul(&x2, irr)).add(b)
}

/// The precomputed meet-in-the-middle pair table: `Pᵢ + Pⱼ` for all
/// `i ≤ j`, keyed by the sum.  Complete by construction, so a lookup
/// miss is a proof that no pair sums to the queried point.
struct PairTable {
    sums: HashMap<FastPoint, (u16, u16)>,
    entries: usize,
    retained_bytes: usize,
    build_adds: u64,
}

impl PairTable {
    fn build(curve: &FastCurve, base: &[FastPoint]) -> Self {
        let mut sums = HashMap::with_capacity(base.len() * base.len() / 2 + 1);
        let mut build_adds = 0u64;
        for i in 0..base.len() {
            for j in i..base.len() {
                let s = curve.add(base[i], base[j]);
                build_adds += 1;
                sums.entry(s).or_insert((i as u16, j as u16));
            }
        }
        let entries = sums.len();
        // FastPoint is Copy { x, y, infinity } (17 bytes padded to 24);
        // the payload adds 4 bytes of indices plus hashmap overhead,
        // accounted per entry as 24 + 8 for the key/value pair.
        let retained_bytes = entries * 32;
        Self {
            sums,
            entries,
            retained_bytes,
            build_adds,
        }
    }
}

/// One target and the two arms' outcomes.
#[derive(Clone, serde::Serialize)]
struct TargetRow {
    kind: &'static str,
    mitm_found: bool,
    pairs_found: bool,
    agree: bool,
    mitm_ms: f64,
    pairs_ms: f64,
    pairs_pairs: u64,
    mitm_ops: u64,
    witness_group_verified: bool,
}

/// Decide, group-exactly, whether `p` or `−p` is a sum of three base
/// points.  Returns the witness indices and the number of point
/// additions performed.
fn mitm_decide(
    curve: &FastCurve,
    base: &[FastPoint],
    table: &PairTable,
    p: FastPoint,
) -> (Option<(usize, usize, usize, bool)>, u64) {
    let mut ops = 0u64;
    for (negated, target) in [false, true]
        .into_iter()
        .map(|neg| (neg, if neg { curve.neg(p) } else { p }))
    {
        for (k, &p3) in base.iter().enumerate() {
            let q = curve.add(target, curve.neg(p3));
            ops += 2;
            if let Some(&(i, j)) = table.sums.get(&q) {
                return (Some((i as usize, j as usize, k, negated)), ops);
            }
        }
    }
    (None, ops)
}

fn main() {
    let (curve, gf, oracle, base, _b_elt) = setup();
    println!("binary MITM oracle: n = {N}, ℓ = {L}, non-Koblitz a = z, b = 1");
    println!(
        "factor base: {} points over the {}-dimensional subspace",
        base.len(),
        L
    );

    // ── factor-base stage: build the pair table ──────────────────────
    let t0 = Instant::now();
    let table = PairTable::build(&curve, &base);
    let construction_wall_ms = t0.elapsed().as_secs_f64() * 1e3;
    println!(
        "pair table: {} distinct sums, {} build adds, {:.1} ms, ~{} KiB",
        table.entries,
        table.build_adds,
        construction_wall_ms,
        table.retained_bytes / 1024
    );

    // ── targets ───────────────────────────────────────────────────────
    let irr = find_irreducible_sparse(N).unwrap();
    let a = F2mElement::from_biguint(&BigUint::from(A_BITS), N);
    let b_elt = F2mElement::from_biguint(&BigUint::from(B), N);
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);

    let mut targets: Vec<(FastPoint, &'static str)> = Vec::new();
    while targets.len() < NATURAL_TARGETS {
        let cand = rng.next() & ((1u64 << N) - 1);
        if cand == 0 {
            continue;
        }
        let x = F2mElement::from_biguint(&BigUint::from(cand), N);
        let c = curve_rhs(&x, &a, &b_elt, &irr);
        let ys = solve_y(&x, &c, &irr);
        if ys.is_empty() {
            continue;
        }
        let y = &ys[(rng.next() & 1) as usize % ys.len()];
        let p = FastPoint::affine(gf.from_element(&x), gf.from_element(y));
        assert!(curve.is_on_curve(p));
        targets.push((p, "natural"));
    }
    for _ in 0..PLANTED_TARGETS {
        // Sum of three distinct base points: satisfiable by construction.
        let i = (rng.next() as usize) % base.len();
        let j = (rng.next() as usize) % base.len();
        let k = (rng.next() as usize) % base.len();
        let p = curve.add(curve.add(base[i], base[j]), base[k]);
        if p.infinity || !curve.is_on_curve(p) {
            continue;
        }
        targets.push((p, "planted_sat"));
    }

    // ── per-target comparison ───────────────────────────────────────
    let mut rows: Vec<TargetRow> = Vec::new();
    let mut disagreements = 0usize;
    for (p, kind) in &targets {
        let (mitm_witness, mitm_ops) = mitm_decide(&curve, &base, &table, *p);
        let mitm_found = mitm_witness.is_some();
        let mut witness_group_verified = false;
        if let Some((i, j, k, negated)) = &mitm_witness {
            let sum = curve.add(curve.add(base[*i], base[*j]), base[*k]);
            let want = if *negated { curve.neg(*p) } else { *p };
            witness_group_verified = sum == want;
            assert!(
                witness_group_verified,
                "MITM witness must re-add to the target"
            );
        }

        // Reference arm, timed on the abscissa of the target point.
        let xr = p.x;
        let t1 = Instant::now();
        let (w, pairs_pairs) = oracle.decompose(xr, &gf);
        let pairs_ms = t1.elapsed().as_secs_f64() * 1e3;
        let pairs_found = w.is_some();

        let t2 = Instant::now();
        let _ = mitm_decide(&curve, &base, &table, *p);
        let mitm_ms = t2.elapsed().as_secs_f64() * 1e3;

        let agree = mitm_found == pairs_found;
        if !agree {
            disagreements += 1;
        }
        rows.push(TargetRow {
            kind,
            mitm_found,
            pairs_found,
            agree,
            mitm_ms,
            pairs_ms,
            pairs_pairs,
            mitm_ops,
            witness_group_verified,
        });
    }

    // ── report ───────────────────────────────────────────────────────
    let mut mitm_times: Vec<f64> = rows.iter().map(|r| r.mitm_ms).collect();
    let mut pairs_times: Vec<f64> = rows.iter().map(|r| r.pairs_ms).collect();
    mitm_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    pairs_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = |v: &[f64]| v[v.len() / 2];
    let found = rows.iter().filter(|r| r.mitm_found).count();
    let refuted = rows.len() - found;
    let ratio = median(&pairs_times) / median(&mitm_times);
    println!(
        "\n{:<8} {:>10} {:>10} {:>7} {:>12} {:>9}",
        "kind", "mitm_ms", "pairs_ms", "agree", "pairs_pairs", "mitm_ops"
    );
    for r in &rows {
        println!(
            "{:<8} {:>10.3} {:>10.3} {:>7} {:>12} {:>9}",
            r.kind, r.mitm_ms, r.pairs_ms, r.agree, r.pairs_pairs, r.mitm_ops
        );
    }
    println!(
        "\nverdict mix: {found} found / {refuted} refuted of {}",
        rows.len()
    );
    println!("disagreements: {disagreements}");
    println!(
        "median MITM: {:.3} ms, median pairs-and-solve: {:.3} ms",
        median(&mitm_times),
        median(&pairs_times)
    );
    println!("speedup (median/median): {ratio:.1}x");
    let gate = median(&mitm_times) <= 0.5 * median(&pairs_times);
    println!(
        "sub-2^(2l) gate (median <= half pairs-and-solve): {}",
        if gate { "PASS" } else { "FAIL" }
    );

    let report = serde_json::json!({
        "regime": "binary",
        "stage": "decomposition",
        "claim": "first sub-2^(2l) decomposition oracle at l=8 on a non-Koblitz binary curve",
        "n": N,
        "l": L,
        "curve": { "model": "y^2+xy=x^3+a*x^2+b", "a": "z", "b": 1, "non_koblitz": true },
        "construction_method": "group-law meet-in-the-middle pair table over subspace factor-base points",
        "oracle_class": "mitm_pair_table",
        "factor_base": {
            "dimension": L,
            "points": base.len(),
            "materialized": true,
            "construction_wall_ms": construction_wall_ms,
            "retained_bytes": table.retained_bytes,
            "table_entries": table.entries,
            "build_adds": table.build_adds,
        },
        "decomposition": {
            "unknowns": format!("3*{} = {} point coordinates (group level)", L, 3 * L),
            "system_degree": "group law (no polynomial system)",
            "oracle_class": "mitm_pair_table",
            "ffd": {
                "status": "unmeasured",
                "reason": "group-law MITM; no polynomial system is built or solved"
            },
            "targets": rows.len(),
            "verdict_mix": { "found": found, "refuted": refuted },
            "median_ms_per_target": { "mitm": median(&mitm_times), "pairs_and_solve": median(&pairs_times) },
            "speedup_median_ratio": ratio,
            "gate_median_le_half_pairs": gate,
            "disagreements": disagreements,
            "per_target": rows,
            "all_witnesses_group_verified": rows.iter().all(|r| !r.mitm_found || r.witness_group_verified),
            "per_target_ops": { "mitm_point_adds_max": rows.iter().map(|r| r.mitm_ops).max() },
        },
        "non_claims": [
            "no key recovery; public synthetic targets only",
            "not an end-to-end ECDLP or vs_rho claim",
            "table build is charged to the factor-base stage, not the per-target oracle",
            "constant-factor win; no asymptotic claim",
        ],
    });
    let out_dir = std::path::Path::new("research/binary_mitm_oracle_l8_20261015");
    std::fs::create_dir_all(out_dir).unwrap();
    let out_path = out_dir.join("results.json");
    std::fs::write(&out_path, serde_json::to_vec_pretty(&report).unwrap()).unwrap();
    println!("\nwrote {}", out_path.display());
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Small-rung control (`n = 15`, `ℓ = 5`): the MITM table's decisions
    /// must agree with an exhaustive enumeration of all base-point triples,
    /// summed in the group, and with the abscissa-level pairs-and-solve
    /// oracle on the same targets.
    #[test]
    fn mitm_agrees_with_exhaustive_triples_and_pairs_and_solve() {
        let irr = find_irreducible_sparse(15).expect("sparse irreducible of degree 15");
        let a = F2mElement::from_biguint(&BigUint::from(A_BITS), 15);
        let b_elt = F2mElement::from_biguint(&BigUint::from(B), 15);
        let mut generator = None;
        let mut probe = 1u64;
        while generator.is_none() {
            let x = F2mElement::from_biguint(&BigUint::from(probe), 15);
            let c = curve_rhs(&x, &a, &b_elt, &irr);
            for y in solve_y(&x, &c, &irr) {
                generator = Some(BinaryPoint::Affine { x, y });
                break;
            }
            probe += 1;
        }
        let curve = BinaryCurve {
            m: 15,
            irreducible: irr.clone(),
            a: a.clone(),
            b: b_elt.clone(),
            generator: generator.unwrap(),
            order: BigUint::one(),
            cofactor: BigUint::one(),
        };
        let fast = FastCurve::new(&curve).unwrap();
        let gf = Gf2::new(&irr);
        let oracle = SubspaceOracle::low_order(5, B, &gf);

        // Factor base over ⟨1, z, …, z⁴⟩.
        let mut base: Vec<FastPoint> = Vec::new();
        for v in 0..(1u64 << 5) {
            let x = F2mElement::from_biguint(&BigUint::from(v), 15);
            let c = curve_rhs(&x, &a, &b_elt, &irr);
            for y in solve_y(&x, &c, &irr) {
                let p = FastPoint::affine(gf.from_element(&x), gf.from_element(&y));
                assert!(fast.is_on_curve(p));
                base.push(p);
            }
        }
        base.sort_by_key(|p| (p.x, p.y));
        base.dedup();
        let table = PairTable::build(&fast, &base);

        // Exhaustive set of all triple sums (as x-level verdict set).
        let mut triple_sums: std::collections::HashSet<FastPoint> = Default::default();
        for i in 0..base.len() {
            for j in i..base.len() {
                let s = fast.add(base[i], base[j]);
                for k in 0..base.len() {
                    triple_sums.insert(fast.add(s, base[k]));
                }
            }
        }

        let mut rng = Rng(0x1234_5678_9ABC_DEF0);
        let mut checked = 0;
        let mut found = 0;
        while checked < 200 {
            let cand = rng.next() & ((1u64 << 15) - 1);
            if cand == 0 {
                continue;
            }
            let x = F2mElement::from_biguint(&BigUint::from(cand), 15);
            let c = curve_rhs(&x, &a, &b_elt, &irr);
            let ys = solve_y(&x, &c, &irr);
            if ys.is_empty() {
                continue;
            }
            let p = FastPoint::affine(gf.from_element(&x), gf.from_element(&ys[0]));
            // Ground truth: the point level (either sign) is a triple sum.
            let exhaustive = triple_sums.contains(&p) || triple_sums.contains(&fast.neg(p));
            // MITM (both signs), with the witness re-added in the group.
            let (w, _ops) = mitm_decide(&fast, &base, &table, p);
            if let Some((i, j, k, negated)) = w {
                let sum = fast.add(fast.add(base[i], base[j]), base[k]);
                let want = if negated { fast.neg(p) } else { p };
                assert_eq!(sum, want, "witness must re-add to the target");
            }
            let mitm = w.is_some();
            assert_eq!(
                exhaustive, mitm,
                "MITM must agree with exhaustive triple enumeration"
            );
            // Pairs-and-solve at the abscissa level.
            let (ps, _pairs) = oracle.decompose(p.x, &gf);
            assert_eq!(exhaustive, ps.is_some(), "pairs-and-solve must agree too");
            if exhaustive {
                found += 1;
            }
            checked += 1;
        }
        assert!(
            found > 0 && found < checked,
            "mix of verdicts expected at n=15, l=5"
        );
    }
}
