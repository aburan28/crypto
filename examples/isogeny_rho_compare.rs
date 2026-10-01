//! Pollard-rho cost across a frozen binary Koblitz isogeny class (CPU screen).
//!
//! Reuses the frozen `a₆` census from
//! `experiments/koblitz_isogeny_cost_sweep.json` (n=17, 273 members). Runs
//! negation-only rho on the Koblitz member and a stratified sample of other
//! members, plus signed-Frobenius rho on the Koblitz member only.
//!
//!     cargo run --release --example isogeny_rho_compare -- \
//!         research/notes/ecc2k130/isogeny_rho_speed_20260930/evidence/n17-rho.json

use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::ic_boundary::{
    binary_point_count, rho_reference_negation, ArtinSchreier, BinaryGroup, CountedGroup, GroupOps,
    RhoResult,
};
use crypto_lib::cryptanalysis::koblitz_fast::{FastCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    find_irreducible_sparse, koblitz_signed_frobenius_rho_reference, KoblitzCurve,
    KoblitzSignedRhoOptions,
};
use crypto_lib::cryptanalysis::koblitz_isogeny_cost::preferred_family;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::fs;
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

fn load_members(n: u32) -> (u8, Vec<u64>, u64) {
    let raw =
        fs::read_to_string("experiments/koblitz_isogeny_cost_sweep.json").expect("frozen census");
    let v: Value = serde_json::from_str(&raw).expect("json");
    for sweep in v["sweeps"].as_array().expect("sweeps") {
        if sweep["n"].as_u64() == Some(n as u64) {
            let a2 = sweep["census"]["a2"].as_u64().unwrap() as u8;
            let order = sweep["census"]["target_order"].as_u64().unwrap();
            let mut members: Vec<u64> = sweep["rows"]
                .as_array()
                .unwrap()
                .iter()
                .map(|r| r["a6"].as_u64().unwrap())
                .collect();
            members.sort_unstable();
            members.dedup();
            return (a2, members, order);
        }
    }
    panic!("no frozen sweep for n={n}");
}

fn sample_members(members: &[u64], k: usize) -> Vec<u64> {
    let mut out = vec![1u64];
    let others: Vec<u64> = members.iter().copied().filter(|&a| a != 1).collect();
    if others.is_empty() {
        return out;
    }
    for i in 1..=k {
        let idx = (i * (others.len().saturating_sub(1))) / k;
        let a = others[idx.min(others.len() - 1)];
        if !out.contains(&a) {
            out.push(a);
        }
    }
    out
}

fn to_element(v: u64, n: u32) -> F2mElement {
    let bits: Vec<u32> = (0..n).filter(|i| (v >> i) & 1 == 1).collect();
    F2mElement::from_bit_positions(&bits, n)
}

fn factorise(mut v: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut d = 2u64;
    while d * d <= v {
        let mut e = 0u32;
        while v.is_multiple_of(d) {
            v /= d;
            e += 1;
        }
        if e > 0 {
            out.push((d, e));
        }
        d += if d == 2 { 1 } else { 2 };
    }
    if v > 1 {
        out.push((v, 1));
    }
    out
}

struct Built {
    fast: FastCurve,
    generator: FastPoint,
    r: u64,
    order: u64,
}

fn rebuild(n: u32, irr: &IrreduciblePoly, a2: u8, a6: u64, seed: u64) -> Built {
    let gf = Gf2::new(irr);
    let ash = ArtinSchreier::new(&gf);
    let order = binary_point_count(&gf, &ash, a2 as u64, a6);
    let &(r, e) = factorise(order).last().expect("factors");
    assert_eq!(e, 1, "repeated largest prime");
    let h = order / r;
    let curve = BinaryCurve {
        m: n,
        irreducible: irr.clone(),
        a: to_element(a2 as u64, n),
        b: to_element(a6, n),
        generator: BinaryPoint::Infinity,
        order: BigUint::from(r),
        cofactor: BigUint::from(h),
    };
    let fast = FastCurve::new(&curve).expect("fast");
    let group = BinaryGroup(&fast);
    let mut ops = GroupOps::default();
    let mut rng = StdRng::seed_from_u64(seed ^ 0xA6A6 ^ a6);
    let mask = if n >= 64 { u64::MAX } else { (1u64 << n) - 1 };
    let generator = loop {
        let x = rng.gen::<u64>() & mask;
        if x == 0 {
            continue;
        }
        let inv = fast.field.inv(x);
        let c = x ^ (a2 as u64) ^ fast.field.mul(a6, fast.field.sqr(inv));
        let Some(u) = ash.solve(c) else {
            continue;
        };
        let y = fast.field.mul(x, u);
        let p = FastPoint::affine(x, y);
        let g = group.mul(&mut ops, p, h);
        if g.infinity {
            continue;
        }
        if !group.mul(&mut ops, g, r).infinity {
            continue;
        }
        break g;
    };
    Built {
        fast,
        generator,
        r,
        order,
    }
}

fn plant_target(
    group: &BinaryGroup<'_>,
    generator: FastPoint,
    r: u64,
    seed: u64,
) -> (u64, FastPoint) {
    let mut rng = StdRng::seed_from_u64(seed);
    let d = rng.gen_range(1..r);
    let mut ops = GroupOps::default();
    let q = group.mul(&mut ops, generator, d);
    (d, q)
}

fn negation_row(n: u32, irr: &IrreduciblePoly, a2: u8, a6: u64, seed: u64, reps: usize) -> Value {
    let me = rebuild(n, irr, a2, a6, seed);
    let group = BinaryGroup(&me.fast);
    let max_steps = ((std::f64::consts::PI * me.r as f64 / 2.0).sqrt() * 200.0) as u64;
    let mut runs = Vec::new();
    for rep in 0..reps {
        let run_seed = seed.wrapping_add(rep as u64 * 0x9E37_79B9);
        let (planted, q) = plant_target(&group, me.generator, me.r, run_seed ^ 0xC0FFEE);
        let report: RhoResult =
            rho_reference_negation(&group, me.generator, q, me.r, run_seed, max_steps);
        runs.push(json!({
            "rep": rep,
            "planted": planted,
            "recovered": report.recovered,
            "verified": report.verified && report.recovered == Some(planted),
            "automorphisms": report.automorphisms,
            "steps": report.steps,
            "gae": report.gae,
            "S": report.s,
            "S_walk": report.s_walk,
            "expected_steps": report.expected_steps,
            "steps_over_expected": report.steps_over_expected,
            "wall_ns": report.wall_ns,
            "method": report.method,
        }));
    }
    json!({
        "a6": a6,
        "is_koblitz": a6 == 1,
        "arm": "negation",
        "order": me.order,
        "r": me.r,
        "A": 2,
        "S_floor": (std::f64::consts::PI / 4.0).sqrt(),
        "runs": runs,
    })
}

fn signed_frobenius_row(n: u32, a2: u8, seed: u64, reps: usize) -> Value {
    let kc = KoblitzCurve::new(a2, n).expect("koblitz");
    let r = kc.subgroup_order.to_u64().expect("r fits u64");
    let g = kc.generator().clone();
    let mut runs = Vec::new();
    for rep in 0..reps {
        let run_seed = seed.wrapping_add(rep as u64 * 0x9E37_79B9);
        let mut rng = StdRng::seed_from_u64(run_seed ^ 0x5F00_EE01);
        let planted = rng.gen_range(1..r);
        let planted_big = BigUint::from(planted);
        let q = kc.mul(&g, &planted_big);
        let opts = KoblitzSignedRhoOptions {
            seed: run_seed,
            ..Default::default()
        };
        let t0 = Instant::now();
        let report = koblitz_signed_frobenius_rho_reference(&kc, &q, &opts, &mut |_| {});
        let wall_ns = t0.elapsed().as_nanos() as u64;
        let recovered = report.recovered_log.as_ref().and_then(|v| v.to_u64());
        let verified = report.verified
            && recovered
                .map(|s| kc.mul(&g, &BigUint::from(s)) == q)
                .unwrap_or(false)
            && recovered == Some(planted);
        let steps = report.charges.walk_group_additions;
        let setup = report.charges.setup_group_additions;
        let bits = (r as f64).log2();
        let smuls = report.charges.setup_scalar_multiplications
            + report.charges.candidate_verification_scalar_multiplications;
        let gae = (setup + steps) as f64 + 1.5 * bits * smuls as f64;
        let s = gae / (r as f64).sqrt();
        let expected = (std::f64::consts::PI * r as f64 / 2.0).sqrt() / (2.0 * n as f64).sqrt();
        runs.push(json!({
            "rep": rep,
            "planted": planted,
            "recovered": recovered,
            "verified": verified,
            "automorphisms": 2 * n,
            "steps": steps,
            "iterations": report.iterations,
            "gae": gae,
            "S": s,
            "S_walk": steps as f64 / (r as f64).sqrt(),
            "expected_steps": expected,
            "steps_over_expected": steps as f64 / expected,
            "wall_ns": wall_ns,
            "method": "signed-Frobenius koblitz_signed_frobenius_rho_reference",
        }));
    }
    json!({
        "a6": 1,
        "is_koblitz": true,
        "arm": "signed_frobenius",
        "r": r,
        "A": 2 * n,
        "S_floor": (std::f64::consts::PI / (4.0 * n as f64)).sqrt(),
        "runs": runs,
    })
}

fn main() {
    let out = PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| {
        "research/notes/ecc2k130/isogeny_rho_speed_20260930/evidence/n17-rho.json".into()
    }));
    let n = 17u32;
    let (a2_pref, _r_pref, _h_pref) = preferred_family(n).expect("family");
    let (a2, members, order) = load_members(n);
    assert_eq!(
        a2, a2_pref,
        "frozen census family must match preferred_family"
    );
    let sample = sample_members(&members, 8);
    let irr = find_irreducible_sparse(n).expect("irr");
    let base_seed = 2026093001u64 + (n as u64) * 1000;
    let reps = 3usize;

    println!(
        "n={n} a2={a2} class_members={} sample={sample:?}",
        members.len()
    );
    let t0 = Instant::now();
    let mut rows = Vec::new();
    println!("signed_frobenius a6=1");
    rows.push(signed_frobenius_row(n, a2, base_seed, reps));
    for &a6 in &sample {
        let seed = base_seed + (a6 % 997);
        println!("negation a6={a6}");
        rows.push(negation_row(n, &irr, a2, a6, seed, reps));
    }

    let payload = json!({
        "schema": "isogeny_rho_compare/v1",
        "protocol": "research/notes/ecc2k130/isogeny_rho_speed_20260930/PROTOCOL.md",
        "generated_unix": SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs(),
        "n": n,
        "a2": a2,
        "class_order": order,
        "class_members": members.len(),
        "sample_a6": sample,
        "reps": reps,
        "host": {
            "rayon_num_threads": std::env::var("RAYON_NUM_THREADS").ok(),
            "seconds": t0.elapsed().as_secs_f64(),
        },
        "rows": rows,
        "prediction": {
            "S_leaf_over_S_E0": (2.0 * n as f64).sqrt(),
            "note": "non-Koblitz members eligible for A=2 only; Koblitz signed-Frobenius A=2n"
        }
    });
    if let Some(parent) = out.parent() {
        fs::create_dir_all(parent).ok();
    }
    fs::write(&out, serde_json::to_string_pretty(&payload).unwrap() + "\n").unwrap();
    println!("wrote {}", out.display());
}
