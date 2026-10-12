//! Measure once-per-curve symbolic setup and repeated S4 specialization separately.
use crypto_lib::cryptanalysis::gaudry_cubic::*;
use serde_json::json;
use std::{env, fs, hint::black_box, time::Instant};

fn main() {
    let a: Vec<String> = env::args().collect();
    assert_eq!(a.len(), 7, "p seed targets setup_repetitions evaluation_repetitions output");
    let p: u64 = a[1].parse().unwrap();
    let seed: u64 = a[2].parse().unwrap();
    let count: u64 = a[3].parse().unwrap();
    let setup_reps: u64 = a[4].parse().unwrap();
    let eval_reps: u64 = a[5].parse().unwrap();
    assert!(count > 0 && setup_reps > 0 && eval_reps > 0);
    let inst = generate_instance3(p, seed);
    let f = &inst.curve.field;
    f.reset_muls();
    let start = Instant::now();
    for _ in 0..setup_reps {
        black_box(SymmetrisedS4::precompute(black_box(&inst.curve)));
    }
    let setup_s = start.elapsed().as_secs_f64();
    let setup_muls = f.muls();
    let pre = SymmetrisedS4::precompute(&inst.curve);
    let targets: Vec<E3> = (0..count).map(|k| E3([
        (k * 7919) % p, (k * 104729 + seed) % p, (k * 15485863 + 7) % p,
    ])).collect();
    // Include zero and base-field targets explicitly, in addition to arbitrary
    // extension-field abscissas (which need not lift to curve points).
    let mut targets = targets;
    targets[0] = Fp3::ZERO;
    if targets.len() > 1 { targets[1] = Fp3::ONE; }
    f.reset_muls();
    let start = Instant::now();
    for _ in 0..eval_reps {
        for target in &targets {
            black_box(pre.weil_restrict(black_box(f), black_box(target)));
        }
    }
    let specialization_s = start.elapsed().as_secs_f64();
    let specialization_muls = f.muls();
    let mut hash = blake3::Hasher::new();
    for target in &targets {
        let polys = pre.weil_restrict(f, target);
        let canonical: Vec<_> = polys.into_iter().map(|poly| {
            let mut terms: Vec<_> = poly.into_iter().collect();
            terms.sort_unstable();
            terms
        }).collect();
        hash.update(&serde_json::to_vec(&(target, canonical)).unwrap());
    }
    fs::write(&a[6], serde_json::to_string_pretty(&json!({
        "p":p,"seed":seed,"targets":count,"setup_repetitions":setup_reps,
        "evaluation_repetitions":eval_reps,"setup_s":setup_s,
        "specialization_s":specialization_s,"setup_fp_muls":setup_muls,
        "specialization_fp_muls":specialization_muls,
        "polynomial_digest":hash.finalize().to_hex().as_str(),
    })).unwrap()).unwrap();
}
