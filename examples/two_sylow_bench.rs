//! Component timings of the 2-Sylow halving search at P-192
//! (`src/cryptanalysis/two_sylow_search.rs`): the primitives each draw
//! uses, per call, median of `--reps` repetitions over 4096 inputs.
//!
//! ```text
//! cargo run --release --example two_sylow_bench -- [--reps 7] [--json OUT]
//! ```

use std::hint::black_box;
use std::time::Instant;

use crypto_lib::cryptanalysis::fp_backend::{FpP192, PrimeField};
use crypto_lib::cryptanalysis::two_sylow_search::{
    draw_rng, element_from_rng, jacobi_p192_fast, p192_pow_isqrt, Arm, Family, Halver,
    HalvingField, SplitMix,
};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .and_then(|i| args.get(i + 1).cloned())
    };
    let reps: usize = get("--reps").and_then(|v| v.parse().ok()).unwrap_or(7);
    let f = FpP192::new();
    let mut rng = SplitMix(42);
    let xs: Vec<[u64; 3]> = (0..4096).map(|_| element_from_rng(&f, &mut rng)).collect();
    let n = xs.len() as f64;
    let mut rows = Vec::new();
    let mut time = |name: &str, body: &mut dyn FnMut() -> u64| {
        let mut v: Vec<f64> = (0..reps)
            .map(|_| {
                let t = Instant::now();
                black_box(body());
                t.elapsed().as_secs_f64() * 1e9 / n
            })
            .collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = v[v.len() / 2];
        println!(
            "{name:<34} median {med:>9.1} ns  (min {:.1}, max {:.1})",
            v[0],
            v[v.len() - 1]
        );
        rows.push(serde_json::json!({"op": name, "median_ns": med, "min_ns": v[0], "max_ns": v[v.len()-1], "reps": reps, "inputs": xs.len()}));
    };
    time("mul", &mut || {
        let mut acc = f.one();
        for &x in &xs {
            acc = f.mul(acc, x);
        }
        acc[0]
    });
    time("sqr (dedicated)", &mut || {
        let mut s = 0u64;
        for &x in &xs {
            s ^= f.sqr(x)[0];
        }
        s
    });
    time("isqrt addition chain", &mut || {
        let mut s = 0u64;
        for &x in &xs {
            s ^= p192_pow_isqrt(&f, x)[0];
        }
        s
    });
    time("sqrt (PrimeField::sqrt, window)", &mut || {
        let mut s = 0u64;
        for &x in &xs {
            s ^= f.sqrt(x).map_or(0, |r| r[0]);
        }
        s
    });
    time("legendre (jacobi_limbs binary)", &mut || {
        let mut s = 0i64;
        for &x in &xs {
            s += i64::from(f.legendre(x));
        }
        s as u64
    });
    time("jacobi posdivsteps (fast)", &mut || {
        let mut s = 0i64;
        for &x in &xs {
            s += i64::from(jacobi_p192_fast(&x));
        }
        s as u64
    });
    time("chi via HalvingField", &mut || {
        let mut s = 0i64;
        for &x in &xs {
            s += i64::from(f.chi(x));
        }
        s as u64
    });
    time("inv (bingcd)", &mut || {
        let mut s = 0u64;
        for &x in &xs {
            s ^= f.inv(x).map_or(0, |r| r[0]);
        }
        s
    });
    time("element_from_rng (draw m)", &mut || {
        let mut s = 0u64;
        let mut r = SplitMix(7);
        for _ in 0..xs.len() {
            s ^= element_from_rng(&f, &mut r)[0];
        }
        s
    });
    let hv = Halver::new(f.clone(), Arm::Fast);
    time("draw_curve z2xz8", &mut || {
        let mut s = 0u64;
        for i in 0..xs.len() as u64 {
            s ^= hv.draw_curve(Family::Z2xZ8, i).map_or(0, |c| c.0.e[1][0]);
        }
        s
    });
    time("draw_rng only", &mut || {
        let mut s = 0u64;
        for i in 0..xs.len() as u64 {
            s ^= draw_rng(Family::Z2xZ8, i).next_u64();
        }
        s
    });
    if let Some(out) = get("--json") {
        std::fs::write(
            out,
            serde_json::to_string_pretty(&serde_json::json!({"components": rows})).unwrap(),
        )
        .unwrap();
    }
}
