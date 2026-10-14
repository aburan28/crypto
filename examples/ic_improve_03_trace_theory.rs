//! IC improvement item 03, theory gate for the trace pre-filter.
//!
//! On `E: y^2 + xy = x^3 + a x^2 + b` over `F_{2^n}` (n odd), the map
//! `chi(P) = Tr(x_P) + Tr(a)` (with `chi(O) = 0`) is the homomorphism
//! `E(F) -> E(F)/2E(F) = Z/2`, and the curve equation gives
//! `chi(P) = Tr(sqrt(b)/x_P)`. Three checks:
//!
//! 1. `identity`: for random `x1, x2` and every root `X` of `S3(x1, x2, X) = 0`,
//!    `Tr(x1) + Tr(x2) = Tr(sqrt(b)/X)` (random `b`, and the Koblitz `b = 1`).
//! 2. `control`: on uniformly drawn curve abscissae (item 08's target law) the
//!    filter `Tr(sqrt(b)/X) = 1` fires on about half; on subgroup points
//!    (`[k]G`, and `[k]G - R` for an orbit representative `R`) it never fires.
//! 3. `exhaustive`: on a small field, for subspace bases `V` (inside `ker Tr`,
//!    inside another hyperplane, random), every m-fold sum (m = 1..4) of points
//!    with abscissa in `V`; tabulates `chi` of the sums and the F2-span
//!    dimension of `{sqrt(b)/X}` and `{X}` over the sums (a linear condition
//!    other than the trace would show as a span below `n - 1`).
//!
//! Usage: `ic_improve_03_trace_theory <out.jsonl>`; prints one JSON line per check.
use crypto_lib::cryptanalysis::koblitz_fast_arith::{s3_x_roots, FastBinaryCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::collections::HashSet;
use std::io::Write;

fn trace(gf: &Gf2, x: u64) -> u32 {
    let mut acc = x;
    let mut p = x;
    for _ in 1..gf.n {
        p = gf.sqr(p);
        acc ^= p;
    }
    debug_assert!(acc <= 1);
    (acc & 1) as u32
}

fn mask(gf: &Gf2) -> u64 {
    (1u64 << gf.n) - 1
}

fn sqrt(gf: &Gf2, x: u64) -> u64 {
    gf.sqr_k(x, gf.n - 1)
}

/// chi via the filter's formula; 0 for the point at infinity. The formula
/// is undefined at the 2-torsion point (0, sqrt b), where chi = Tr(a)
/// (callers pass Tr(a) as `tr_a`).
fn chi_inv(gf: &Gf2, sqrt_b: u64, tr_a: u32, p: FastPoint) -> u32 {
    match p {
        None => 0,
        Some((0, _)) => tr_a,
        Some((x, _)) => trace(gf, gf.mul(sqrt_b, gf.inv(x))),
    }
}

fn chi_direct(gf: &Gf2, a: u64, p: FastPoint) -> u32 {
    match p {
        None => 0,
        Some((x, _)) => trace(gf, x) ^ trace(gf, a),
    }
}

/// F2 rank of a set of words.
fn span_dim(words: impl Iterator<Item = u64>) -> usize {
    let mut basis = [0u64; 64];
    let mut dim = 0;
    for mut w in words {
        for bit in (0..64).rev() {
            if (w >> bit) & 1 == 0 {
                continue;
            }
            if basis[bit] == 0 {
                basis[bit] = w;
                dim += 1;
                break;
            }
            w ^= basis[bit];
        }
    }
    dim
}

fn emit(out: &mut std::fs::File, v: Value) {
    println!("{v}");
    writeln!(out, "{v}").unwrap();
}

fn check_identity(out: &mut std::fs::File, n: u32, samples: u64, rng: &mut StdRng) {
    let curve = KoblitzCurve::new(0, n).unwrap();
    let gf = Gf2::new(&curve.curve.irreducible);
    for b_kind in ["koblitz_b1", "random_b"] {
        let (mut roots_checked, mut violations, mut no_root) = (0u64, 0u64, 0u64);
        for _ in 0..samples {
            let b = if b_kind == "koblitz_b1" {
                1
            } else {
                loop {
                    let b = rng.gen::<u64>() & mask(&gf);
                    if b != 0 {
                        break b;
                    }
                }
            };
            let x1 = rng.gen::<u64>() & mask(&gf);
            let x2 = rng.gen::<u64>() & mask(&gf);
            if x1 == 0 || x2 == 0 || x1 == x2 {
                continue;
            }
            let Some(roots) = s3_x_roots(&gf, b, x1, x2) else {
                no_root += 1;
                continue;
            };
            let sb = sqrt(&gf, b);
            for x in roots {
                if x == 0 {
                    continue;
                }
                roots_checked += 1;
                let lhs = trace(&gf, x1) ^ trace(&gf, x2);
                let rhs = trace(&gf, gf.mul(sb, gf.inv(x)));
                if lhs != rhs {
                    violations += 1;
                }
            }
        }
        emit(
            out,
            json!({"check":"identity", "n":n, "b":b_kind, "samples":samples,
                   "pairs_without_root":no_root, "roots_checked":roots_checked,
                   "violations":violations}),
        );
    }
}

fn check_control(out: &mut std::fs::File, n: u32, a: u8, samples: u64, rng: &mut StdRng) {
    let Some(curve) = KoblitzCurve::new(a, n) else {
        emit(out, json!({"check":"control_skipped", "n":n, "a":a, "reason":"no admitted Koblitz rung"}));
        return;
    };
    let fast = FastBinaryCurve::new(&curve.curve.irreducible, a as u64).unwrap();
    let gf = &fast.gf;
    let b = gf.from_element(&curve.curve.b);
    let sb = sqrt(gf, b);
    let tr_a = trace(gf, fast.a);
    let r = curve.subgroup_order.to_u64().unwrap();
    // Uniform curve abscissae.
    let (mut drawn, mut fires, mut chi_mismatch) = (0u64, 0u64, 0u64);
    while drawn < samples {
        let x = rng.gen::<u64>() & mask(gf);
        if x == 0 {
            continue;
        }
        let pts = fast.points_with_x(b, x);
        if pts.is_empty() {
            continue;
        }
        drawn += 1;
        let c = chi_inv(gf, sb, tr_a, pts[0]);
        fires += c as u64;
        if c != chi_direct(gf, fast.a, pts[0]) {
            chi_mismatch += 1;
        }
    }
    emit(
        out,
        json!({"check":"control_uniform_abscissae", "n":n, "a":a, "samples":drawn,
               "filter_fires":fires, "fraction":fires as f64 / drawn as f64,
               "chi_formula_mismatch":chi_mismatch}),
    );
    // Subgroup points [k]G and [k]G - R for a cofactor-projected R.
    let g: FastPoint = match curve.generator() {
        crypto_lib::binary_ecc::BinaryPoint::Affine { x, y } => Some((fast.word(&x), fast.word(&y))),
        _ => unreachable!(),
    };
    let rep = (2u64..)
        .find_map(|x| {
            fast.points_with_x(b, x)
                .into_iter()
                .find_map(|p| fast.scalar_mul(p, &curve.cofactor))
        })
        .map(Some)
        .unwrap();
    let sub_samples = samples / 50;
    let (mut fires_g, mut fires_gr) = (0u64, 0u64);
    for _ in 0..sub_samples {
        let k = rng.gen::<u64>() % (r - 1) + 1;
        let p = fast.scalar_mul(g, &BigUint::from(k));
        fires_g += chi_inv(gf, sb, tr_a, p) as u64;
        let q = fast.add(p, FastBinaryCurve::neg(rep));
        fires_gr += chi_inv(gf, sb, tr_a, q) as u64;
    }
    emit(
        out,
        json!({"check":"control_subgroup", "n":n, "a":a, "samples":sub_samples,
               "filter_fires_kG":fires_g, "filter_fires_kG_minus_R":fires_gr}),
    );
}

fn check_exhaustive(out: &mut std::fs::File, n: u32, a: u8, k: usize, rng: &mut StdRng) {
    let Some(curve) = KoblitzCurve::new(a, n) else {
        emit(out, json!({"check":"exhaustive_skipped", "n":n, "a":a, "reason":"no admitted Koblitz rung"}));
        return;
    };
    let fast = FastBinaryCurve::new(&curve.curve.irreducible, a as u64).unwrap();
    let gf = &fast.gf;
    let b = gf.from_element(&curve.curve.b);
    let sb = sqrt(gf, b);
    let tr_a = trace(gf, fast.a);
    // Whole-group census.
    let mut all = Vec::new();
    for x in 0..=mask(gf) {
        all.extend(fast.points_with_x(b, x));
    }
    let order = all.len() as u64 + 1;
    let all_chi1 = all.iter().filter(|&&p| chi_inv(gf, sb, tr_a, p) == 1).count();
    let sub_chi1 = all
        .iter()
        .filter(|&&p| fast.scalar_mul(p, &curve.subgroup_order).is_none())
        .filter(|&&p| chi_inv(gf, sb, tr_a, p) == 1)
        .count();
    emit(
        out,
        json!({"check":"exhaustive_group", "n":n, "a":a, "group_order":order,
               "points_chi1":all_chi1, "subgroup_points_chi1":sub_chi1,
               "subgroup_order":curve.subgroup_order.to_u64()}),
    );
    // Hyperplanes: ker Tr, and {Tr(c x) = 0} for a random c not in F2.
    let c_other = loop {
        let c = rng.gen::<u64>() & mask(gf);
        if c > 1 {
            break c;
        }
    };
    let in_hyper = |c: u64, x: u64| trace(gf, gf.mul(c, x)) == 0;
    let kinds: [(&str, Option<u64>); 3] = [
        ("inside_ker_Tr", Some(1)),
        ("inside_other_hyperplane", Some(c_other)),
        ("random", None),
    ];
    for (kind, hyper) in kinds {
        for trial in 0..3 {
            // Random k-dim subspace, optionally inside the hyperplane.
            let mut basis: Vec<u64> = Vec::new();
            while basis.len() < k {
                let v = rng.gen::<u64>() & mask(gf);
                if v == 0 || hyper.is_some_and(|c| !in_hyper(c, v)) {
                    continue;
                }
                if span_dim(basis.iter().copied().chain([v])) == basis.len() + 1 {
                    basis.push(v);
                }
            }
            let mut base: Vec<FastPoint> = Vec::new();
            for mask in 1u64..(1 << k) {
                let x = (0..k).filter(|&i| (mask >> i) & 1 == 1).fold(0, |acc, i| acc ^ basis[i]);
                base.extend(fast.points_with_x(b, x));
            }
            let base_chi: Vec<u32> = base.iter().map(|&p| chi_inv(gf, sb, tr_a, p)).collect();
            let base_chi1 = base_chi.iter().filter(|&&c| c == 1).count();
            let mut level: HashSet<FastPoint> = base.iter().copied().collect();
            let mut rows = Vec::new();
            for m in 1..=4usize {
                if m > 1 {
                    let mut next = HashSet::with_capacity(level.len() * 2);
                    for &p in &level {
                        for &q in &base {
                            next.insert(fast.add(p, q));
                        }
                    }
                    level = next;
                }
                let affine: Vec<u64> = level.iter().filter_map(|p| p.map(|(x, _)| x)).collect();
                let chi1 = level.iter().filter(|&&p| chi_inv(gf, sb, tr_a, p) == 1).count();
                let inv_span = span_dim(
                    affine.iter().filter(|&&x| x != 0).map(|&x| gf.mul(sb, gf.inv(x))),
                );
                let x_span = span_dim(affine.iter().copied());
                // Predicted: if every base point has chi = c, every m-sum has chi = m c.
                let uniform = base_chi1 == 0 || base_chi1 == base.len();
                let predicted = uniform.then(|| (m as u32 * base_chi[0]) & 1);
                rows.push(json!({"m":m, "sums":level.len(), "sums_chi1":chi1,
                    "predicted_chi_all":predicted,
                    "prediction_holds":predicted.map(|c| if c == 0 { chi1 == 0 } else { chi1 == level.len() }),
                    "span_dim_sqrtb_over_x":inv_span, "span_dim_x":x_span}));
            }
            emit(
                out,
                json!({"check":"exhaustive_subspace", "n":n, "a":a, "k":k, "kind":kind,
                       "trial":trial, "base_points":base.len(), "base_points_chi1":base_chi1,
                       "by_m":rows}),
            );
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 2, "usage: ic_improve_03_trace_theory <out.jsonl>");
    let mut out = std::fs::File::create(&args[1]).unwrap();
    let mut rng = StdRng::seed_from_u64(0x1C03);
    check_identity(&mut out, 41, 1_000_000, &mut rng);
    for (n, a) in [(41u32, 0u8), (41, 1), (43, 1), (47, 1), (53, 1), (61, 0)] {
        check_control(&mut out, n, a, 1_000_000, &mut rng);
    }
    for (n, a, k) in [(11u32, 0u8, 5usize), (11, 1, 5), (13, 0, 5), (13, 1, 5), (13, 0, 6), (17, 1, 6)] {
        check_exhaustive(&mut out, n, a, k, &mut rng);
    }
}
