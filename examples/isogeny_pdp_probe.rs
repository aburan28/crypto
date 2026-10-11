//! Does an isogenous curve give an easier point-decomposition problem?
//!
//! Toy-size measurement over binary fields, on public synthetic curves.
//! For the Koblitz curve `K_0: y² + xy = x³ + 1` at `n ∈ {17, 19, 23}` and
//! each of its `ℓ`-isogenous neighbours (`ℓ = 3, 5`, found through the
//! modular polynomial `Φ_ℓ mod 2` in `cryptanalysis::binary_isogeny`),
//! with the subspace factor base `V = ⟨1, z, …, z^{l−1}⟩`, it records:
//!
//! 1. `|F|` — how many `x ∈ V` lie on the curve (the on-curve fraction is
//!    `b`-dependent, so this can differ between isogenous curves);
//! 2. exact 2-sum **relation yield**: the fraction of `x ∈ F_{2^n}` that are
//!    `x(P₁ ± P₂)` for `P₁, P₂ ∈ F`, by enumerating every pair;
//! 3. the **Macaulay rank profile and first-fall degree** of the
//!    Weil-descended `S₃` system (`ffd_harness::measure_one`): since `b`
//!    enters the binary `S₃` only as a constant, the prediction is that the
//!    ranks are identical across the isogeny class;
//! 4. the **F4 decomposition cost** (`koblitz_groebner`) on 2-decomposable
//!    targets and on random targets, with the same subspace;
//! 5. whether the model keeps the Frobenius symmetry (`b ∈ F_2`), which is
//!    what gives Koblitz factor bases their `τ`-orbit multiplier.
//!
//! ```bash
//! cargo run --release --example isogeny_pdp_probe -- [--n 17,19,23] [--ell 3,5] [--json out.json]
//! ```

use std::collections::HashSet;
use std::time::Instant;

use crypto_lib::binary_ecc::curve::{point_add, BinaryCurve, BinaryPoint};
use crypto_lib::binary_ecc::f2m::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_isogeny::l_isogenous_neighbours;
use crypto_lib::cryptanalysis::ffd_harness::measure_one;
use crypto_lib::cryptanalysis::ghs_descent::ECurve;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, solve_boolean_system, FieldStructure, SolveOptions,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{find_irreducible, points_with_x};
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;

fn bits_of(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

fn elem_from_bits(bits: u64, n: u32) -> F2mElement {
    let positions: Vec<u32> = (0..n).filter(|k| (bits >> k) & 1 == 1).collect();
    F2mElement::from_bit_positions(&positions, n)
}

/// Arithmetic-only curve object (generator/order are not used here).
fn arith_curve(n: u32, irr: &IrreduciblePoly, a: &F2mElement, b: &F2mElement) -> BinaryCurve {
    BinaryCurve {
        m: n,
        irreducible: irr.clone(),
        a: a.clone(),
        b: b.clone(),
        generator: BinaryPoint::Infinity,
        order: BigUint::from(1u8),
        cofactor: BigUint::from(1u8),
    }
}

struct CurveRow {
    label: String,
    b_bits: u64,
    a_bits: u64,
    frobenius_symmetric: bool,
    fb_size: usize,
    yield2: f64,
    fall_degree: Option<u32>,
    rank_profile: Vec<(u32, u64, u64)>, // (degree, rank, rank_generic)
    f4_ms_decomposable: f64,
    f4_ms_random: f64,
    f4_reductions: f64,
    f4_found: usize,
    f4_targets: usize,
}

fn measure_curve(label: &str, n: u32, irr: &IrreduciblePoly, a: &F2mElement, b: &F2mElement, l: u32, seed: u64) -> CurveRow {
    let curve = arith_curve(n, irr, a, b);
    // Factor base: x in V = span(1, z, ..., z^{l-1}) that lie on the curve.
    let basis: Vec<F2mElement> = (0..l).map(|k| F2mElement::from_bit_positions(&[k], n)).collect();
    let mut fb: Vec<BinaryPoint> = Vec::new();
    for v in 0..(1u64 << l) {
        let x = elem_from_bits(v, n);
        for pt in points_with_x(&curve, &x) {
            fb.push(pt);
        }
    }
    // Exact 2-sum yield: every x(P1 + P2) over unordered pairs (negatives are included since fb holds both points over each x).
    let mut sums: HashSet<u64> = HashSet::new();
    for i in 0..fb.len() {
        for j in i..fb.len() {
            if let BinaryPoint::Affine { x, .. } = point_add(&curve, &fb[i], &fb[j]) {
                sums.insert(bits_of(&x));
            }
        }
    }
    let yield2 = sums.len() as f64 / 2f64.powi(n as i32);
    // Macaulay profile of the full Weil-descended S3 for one random x3.
    let mut rng = StdRng::seed_from_u64(seed);
    let x3 = elem_from_bits(rng.gen_range(1..(1u64 << n)), n);
    let ffd = measure_one(n, irr, b, &x3, 4);
    let rank_profile: Vec<(u32, u64, u64)> = ffd.per_degree.iter().map(|m| (m.degree, m.rank, m.rank_generic)).collect();
    // F4 decomposition cost with the subspace factor base, m = 2.
    let st = FieldStructure::new(n, irr);
    let opts = SolveOptions::default();
    let mut t_dec = 0.0;
    let mut t_rnd = 0.0;
    let mut reductions = 0usize;
    let mut found = 0usize;
    let targets = 6usize;
    for t in 0..targets {
        // decomposable target: R = P_i + P_j
        let i = rng.gen_range(0..fb.len());
        let j = rng.gen_range(0..fb.len());
        let r = point_add(&curve, &fb[i], &fb[j]);
        let BinaryPoint::Affine { x: xr, .. } = r else { continue };
        if let Some(sys) = build_decomposition_system(&basis, &xr, b, 2, &st) {
            let t0 = Instant::now();
            let (sols, stats) = solve_boolean_system(&sys.equations, sys.n_vars, &opts);
            t_dec += t0.elapsed().as_secs_f64() * 1e3;
            reductions += stats.reductions;
            if !sols.is_empty() {
                found += 1;
            }
        }
        // random target
        let xr2 = elem_from_bits(rng.gen_range(1..(1u64 << n)) ^ (t as u64), n);
        if let Some(sys) = build_decomposition_system(&basis, &xr2, b, 2, &st) {
            let t0 = Instant::now();
            let (_sols, stats) = solve_boolean_system(&sys.equations, sys.n_vars, &opts);
            t_rnd += t0.elapsed().as_secs_f64() * 1e3;
            reductions += stats.reductions;
        }
    }
    CurveRow {
        label: label.to_string(),
        b_bits: bits_of(b),
        a_bits: bits_of(a),
        frobenius_symmetric: bits_of(b) <= 1 && bits_of(a) <= 1,
        fb_size: fb.len(),
        yield2,
        fall_degree: ffd.fall_degree,
        rank_profile,
        f4_ms_decomposable: t_dec / targets as f64,
        f4_ms_random: t_rnd / targets as f64,
        f4_reductions: reductions as f64 / (2 * targets) as f64,
        f4_found: found,
        f4_targets: targets,
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let mut ns = vec![17u32, 19, 23];
    let mut ells = vec![3u32, 5];
    let mut json_path: Option<String> = None;
    let mut base_seed: Option<u64> = None; // --base random:<seed> uses a random-b, a=0 base curve instead of K_0
    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "--n" => {
                i += 1;
                ns = argv[i].split(',').filter_map(|s| s.parse().ok()).collect();
            }
            "--ell" => {
                i += 1;
                ells = argv[i].split(',').filter_map(|s| s.parse().ok()).collect();
            }
            "--json" => {
                i += 1;
                json_path = Some(argv[i].clone());
            }
            "--base" => {
                i += 1;
                base_seed = argv[i].strip_prefix("random:").and_then(|v| v.parse().ok());
            }
            _ => {}
        }
        i += 1;
    }
    let mut out = Vec::new();
    println!("# Isogenous curves and the point-decomposition problem (binary toy sizes)\n");
    for &n in &ns {
        let Some(irr) = find_irreducible(n) else { continue };
        let l = n / 2;
        let one = F2mElement::one(n);
        let zero = F2mElement::zero(n);
        let base_b = match base_seed {
            Some(sd) => {
                let mut r = StdRng::seed_from_u64(sd ^ n as u64);
                elem_from_bits(r.gen_range(2..(1u64 << n)), n)
            }
            None => one.clone(),
        };
        let base_label = if base_seed.is_some() { "random-b base" } else { "K_0 (base)" };
        let base = ECurve::new(n, irr.clone(), zero.clone(), base_b.clone());
        let mut curves: Vec<(String, F2mElement, F2mElement)> = vec![(base_label.into(), zero.clone(), base_b.clone())];
        for &ell in &ells {
            for (k, nb) in l_isogenous_neighbours(&base, ell).into_iter().enumerate() {
                // `l_isogenous_neighbours` returns both a-choices per neighbouring j; only the
                // one with the same trace as the base is isogenous, the other is its twist.
                let same_a = nb.a == base.a;
                let kind = if bits_of(&nb.b) == bits_of(&base.b) && same_a { "self (dual isogeny)" } else if same_a { "isogenous" } else { "twist of neighbour" };
                curves.push((format!("{ell}-{kind} #{k}"), nb.a.clone(), nb.b.clone()));
            }
        }
        println!("## n = {n}, factor base V = span(1..z^{}) (2^{l} x-values), 2-decomposition\n", l - 1);
        println!("| curve | a | b | Frobenius-symmetric | \\|F\\| | 2-sum yield | FFD | Macaulay rank (deg 2/3/4) vs generic | F4 ms decomposable | F4 ms random | F4 found | mean reductions |");
        println!("|---|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|");
        for (k, (label, a, b)) in curves.iter().enumerate() {
            let row = measure_curve(label, n, &irr, a, b, l, 1000 + k as u64);
            let prof: Vec<String> = row.rank_profile.iter().map(|(d, r, g)| format!("d{d}: {r}/{g}")).collect();
            println!(
                "| {} | {:#x} | {:#x} | {} | {} | {:.4} | {} | {} | {:.1} | {:.1} | {}/{} | {:.1} |",
                row.label,
                row.a_bits,
                row.b_bits,
                row.frobenius_symmetric,
                row.fb_size,
                row.yield2,
                row.fall_degree.map(|d| d.to_string()).unwrap_or("—".into()),
                prof.join(", "),
                row.f4_ms_decomposable,
                row.f4_ms_random,
                row.f4_found,
                row.f4_targets,
                row.f4_reductions
            );
            out.push(json!({"n":n,"l":l,"curve":row.label,"a":row.a_bits,"b":row.b_bits,"frobenius_symmetric":row.frobenius_symmetric,
                "fb_size":row.fb_size,"yield2":row.yield2,"fall_degree":row.fall_degree,"rank_profile":row.rank_profile,
                "f4_ms_decomposable":row.f4_ms_decomposable,"f4_ms_random":row.f4_ms_random,"f4_found":row.f4_found,"f4_targets":row.f4_targets,
                "f4_mean_reductions":row.f4_reductions}));
        }
        println!();
    }
    if let Some(p) = json_path {
        std::fs::write(&p, serde_json::to_string_pretty(&json!({"results": out})).unwrap()).unwrap();
        println!("JSON written to {p}");
    }
}
