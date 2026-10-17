//! Build the PKM isogeny tower `V = x(R' + ⟨Q⟩)` at the P-192 prime on an
//! auxiliary curve found by `pkm-aux-search`, and check it.
//!
//! ```text
//! cargo run --release --example p192_isogeny_tower -- ARTIFACT.json T [--seed S] [--json OUT]
//! ```
//!
//! The artifact is `backlog_check_aux2t`'s format with `t ≥ T`; the point
//! used is `Q_T = [2^(t−T)]Q`, of order `2^T`. Checks: `|V| = 2^T`
//! distinct abscissae; `L(x) = c` on every `x ∈ V`, with the field
//! backend and, on every `x`, again over `BigUint`; `L(x) ≠ c` on 4096
//! random abscissae outside `V`. Reports how many `x ∈ V` lift to the
//! NIST P-192 curve (`x³ − 3x + b₁₉₂` a square), i.e. the base's size.

use std::time::Instant;

use crypto_lib::cryptanalysis::fp_backend::{FpP192, PrimeField};
use crypto_lib::cryptanalysis::two_sylow_search::{
    big_double, build_tower, element_from_rng, p192, tower_eval, SplitMix, Sw,
};
use num_bigint::BigUint;
use serde_json::{json, Value};

fn num(v: &Value, k: &str) -> BigUint {
    BigUint::parse_bytes(v[k].as_str().expect(k).as_bytes(), 10).expect(k)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let art: Value = serde_json::from_str(&std::fs::read_to_string(&args[1]).unwrap()).unwrap();
    let t: u32 = args[2].parse().unwrap();
    let get = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .and_then(|i| args.get(i + 1).cloned())
    };
    let seed: u64 = get("--seed").and_then(|v| v.parse().ok()).unwrap_or(1);
    let p = p192();
    assert_eq!(num(&art, "p"), p);
    let (a, b) = (num(&art, "a"), num(&art, "b"));
    let ta = art["t"].as_u64().unwrap() as u32;
    assert!(t <= ta && t <= 24);
    let mut q = Some((num(&art, "Qx"), num(&art, "Qy")));
    for _ in 0..ta - t {
        q = big_double(&p, &a, &q);
    }
    let (qx, qy) = q.expect("Q_T");
    let f = FpP192::new();
    let sw = Sw {
        f: f.clone(),
        a: f.from_biguint(&a),
        b: f.from_biguint(&b),
    };
    let start = Instant::now();
    let tower =
        build_tower(&sw, (f.from_biguint(&qx), f.from_biguint(&qy)), t, seed).expect("tower build");
    let build_s = start.elapsed().as_secs_f64();
    let mut v = tower.v.clone();
    v.sort();
    v.dedup();
    let distinct = v.len();
    let fp_ok = tower
        .v
        .iter()
        .all(|&x| tower_eval(&f, &tower.steps, x) == Some(tower.c));
    // BigUint re-evaluation of the composed Vélu x-map.
    let c = f.to_biguint(tower.c);
    let steps: Vec<(BigUint, BigUint)> = tower
        .steps
        .iter()
        .map(|s| (f.to_biguint(s.xi), f.to_biguint(s.tau)))
        .collect();
    let big_ok = tower.v.iter().all(|&x| {
        let mut y = f.to_biguint(x);
        for (xi, tau) in &steps {
            let d = (&y + &p - xi) % &p;
            y = (&y + tau * d.modpow(&(&p - 2u32), &p)) % &p;
        }
        y == c
    });
    let mut rng = SplitMix(seed ^ 0xABCD);
    let set: std::collections::HashSet<[u64; 3]> = tower.v.iter().copied().collect();
    let mut outside = 0;
    let mut outside_hits = 0;
    while outside < 4096 {
        let x = element_from_rng(&f, &mut rng);
        if set.contains(&x) {
            continue;
        }
        outside += 1;
        if tower_eval(&f, &tower.steps, x) == Some(tower.c) {
            outside_hits += 1;
        }
    }
    // The P-192 curve y² = x³ − 3x + b.
    let b192 =
        BigUint::parse_bytes(b"64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1", 16).unwrap();
    let (m3, bb) = (f.neg(f.from_u64(3)), f.from_biguint(&b192));
    let lifting = tower
        .v
        .iter()
        .filter(|&&x| f.legendre(f.add(f.mul(f.add(f.sqr(x), m3), x), bb)) >= 0)
        .count();
    let out = json!({
        "artifact": args[1], "t": t, "artifact_t": ta, "seed": seed,
        "v_len": tower.v.len(), "v_distinct": distinct, "eval_all_c_fp": fp_ok,
        "eval_all_c_biguint": big_ok, "outside_tested": outside, "outside_hits": outside_hits,
        "v_lifting_to_p192": lifting, "build_seconds": build_s,
        "c": c.to_string(),
    });
    println!("{out}");
    if let Some(o) = get("--json") {
        std::fs::write(o, serde_json::to_string_pretty(&out).unwrap() + "\n").unwrap();
    }
    assert!(distinct == 1 << t && fp_ok && big_ok && outside_hits == 0);
}
