//! # Fail-closed validation of the GPU root-finding decomposition oracle
//!
//! Builds an `a = −3` prime-order curve at `BITS` bits
//! (`prime_fast::find_a3_curve`), a small-x factor base of size `FB`,
//! and a target batch of
//!
//! - **positive controls** `R = ±F_i ± F_j` (random `i, j`, random signs),
//!   which the oracle must decompose;
//! - **degenerate controls** `R = ±F_i`;
//! - **random on-curve points** (the rejection path that dominates in
//!   practice);
//! - **production targets**: the `R = aG + bQ` points behind the
//!   relations `prime_fast::collect_relations` actually emitted on this
//!   curve, whose 1- and 2-entry rows the oracle must reproduce;
//!
//! then runs the oracle through `gpu_rootfind::decompose64` (this crate,
//! `gpu/rootfind/rs/`) with every
//! (lookup, chunk) configuration and asserts, for each:
//!
//! 1. every returned relation re-adds to its target in the group;
//! 2. the returned relation *set* equals the CPU `prime_fast` reference
//!    (`crypto_lib::cryptanalysis::gpu_rootfind::cpu_reference64`) on
//!    every target — not just the
//!    controls, the random rejections too (no spurious hits);
//! 3. every positive, degenerate and production control is present.
//!
//! Any failure exits non-zero.  On this machine the oracle runs through
//! the host-emulation backend (the kernel's per-thread function executed
//! on the CPU); with `GPU_ROOTFIND_CUDA_LIB_DIR` set at build time and a
//! device present, `BACKEND=cuda` runs the same checks against the GPU.
//!
//! ```bash
//! # from the repository root
//! cargo run --release --manifest-path gpu/rootfind/rs/Cargo.toml --bin gpu_rootfind_validate
//! BITS=61 FB=4096 RAND_TARGETS=4096 BENCH_OUT=research/gpu_rootfind_decomposition_20261012 \
//!   cargo run --release --manifest-path gpu/rootfind/rs/Cargo.toml --bin gpu_rootfind_validate
//! # knobs: BITS (40..=61) FB RAND_TARGETS CONTROLS BACKEND=auto|emulate|cuda SEED
//! ```
//!
//! ## Not a security result
//!
//! Only the per-target decomposition subroutine is exercised, on public
//! synthetic known-answer instances.  Relation collection and linear
//! algebra are untouched and the asymptotics of index calculus are
//! unchanged; no key-recovery, sub-rho or deployed-curve claim is made.

use std::collections::BTreeSet;
use std::time::Instant;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crypto_lib::cryptanalysis::prime_fast::{collect_relations, find_a3_curve, FastFactorBase, Pt};
use gpu_rootfind::{
    backend_info, cpu_reference64, cuda_available, decompose64, factor_base_canonical,
    targets_canonical, verify_relation64, Backend, Lookup, Options, Relation,
};

fn env_usize(k: &str, d: usize) -> usize {
    std::env::var(k)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(d)
}
fn env_u64(k: &str, d: u64) -> u64 {
    std::env::var(k)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(d)
}

struct Expect {
    target: usize,
    rel: Relation,
    kind: &'static str,
}

fn main() {
    let bits = env_usize("BITS", 48) as u32;
    let fb_size = env_usize("FB", 1024);
    let n_random = env_usize("RAND_TARGETS", 2048);
    let n_controls = env_usize("CONTROLS", 256);
    let seed = env_u64("SEED", 20261012);
    let backend = match std::env::var("BACKEND").as_deref() {
        Ok("cuda") => Backend::Cuda,
        Ok("emulate") => Backend::Emulate,
        _ => Backend::Auto,
    };
    assert!((40..=61).contains(&bits), "BITS must be in 40..=61");

    println!(
        "gpu_rootfind_validate: backend library = {}",
        backend_info()
    );
    println!("  cuda_available = {}", cuda_available());

    let t0 = Instant::now();
    let (curve, g) = find_a3_curve(bits, 3).expect("a=-3 curve");
    println!(
        "  curve {} : p = {} ({} bits), n = {}  [{:.2}s]",
        curve.name,
        curve.f.p,
        bits,
        curve.n,
        t0.elapsed().as_secs_f64()
    );
    let fb = FastFactorBase::build(&curve, fb_size);
    println!("  factor base |F| = {} (max x = {})", fb.len(), fb.max_x);

    // ── targets ──────────────────────────────────────────────────
    let mut rng = StdRng::seed_from_u64(seed);
    let mut targets: Vec<Pt> = Vec::new();
    let mut expect: Vec<Expect> = Vec::new();

    while expect.len() < n_controls {
        let i = rng.gen_range(0..fb.len());
        let j = rng.gen_range(0..fb.len());
        let ei: i8 = if rng.gen_bool(0.5) { 1 } else { -1 };
        let ej: i8 = if rng.gen_bool(0.5) { 1 } else { -1 };
        let pi = if ei < 0 {
            curve.neg(fb.pts[i])
        } else {
            fb.pts[i]
        };
        let pj = if ej < 0 {
            curve.neg(fb.pts[j])
        } else {
            fb.pts[j]
        };
        let Some(r) = curve.add(Some(pi), Some(pj)) else {
            continue;
        };
        let (i, j, ei, ej) = if i <= j {
            (i, j, ei, ej)
        } else {
            (j, i, ej, ei)
        };
        targets.push(r);
        expect.push(Expect {
            target: targets.len() - 1,
            rel: Relation {
                target: (targets.len() - 1) as u32,
                i: i as u32,
                j: j as u32,
                eps: ei,
                tau: ej,
            },
            kind: "positive",
        });
    }
    for _ in 0..16 {
        let i = rng.gen_range(0..fb.len());
        let e: i8 = if rng.gen_bool(0.5) { 1 } else { -1 };
        targets.push(if e < 0 {
            curve.neg(fb.pts[i])
        } else {
            fb.pts[i]
        });
        expect.push(Expect {
            target: targets.len() - 1,
            rel: Relation {
                target: (targets.len() - 1) as u32,
                i: i as u32,
                j: i as u32,
                eps: e,
                tau: 0,
            },
            kind: "degenerate",
        });
    }
    let n_ctl = targets.len();
    while targets.len() < n_ctl + n_random {
        let x = rng.gen_range(1..curve.f.p);
        if let Some(p) = curve.lift_x(x) {
            targets.push(p);
        }
    }
    // production targets: what collect_relations found on this curve
    let q = curve.mul(Some(g), rng.gen_range(1..curve.n)).expect("Q");
    let (rels, stats) = collect_relations(&curve, g, q, &fb, None, 24, 200_000, seed ^ 0xABCD);
    let mut n_prod = 0usize;
    for rel in &rels {
        if rel.entries.len() > 2 {
            continue;
        }
        let r = curve.combine(g, q, rel.a, rel.b).expect("R = aG + bQ");
        targets.push(r);
        let k = targets.len() - 1;
        let want = match rel.entries.as_slice() {
            [(i, s)] => {
                // ±F_i, or 2·F_i (s = ±2) which the oracle reports as (i, s', i, s')
                if s.unsigned_abs() == 2 {
                    let sg = if *s < 0 { -1 } else { 1 };
                    Relation {
                        target: k as u32,
                        i: *i,
                        j: *i,
                        eps: sg,
                        tau: sg,
                    }
                } else {
                    Relation {
                        target: k as u32,
                        i: *i,
                        j: *i,
                        eps: *s,
                        tau: 0,
                    }
                }
            }
            [(i, si), (j, sj)] => Relation {
                target: k as u32,
                i: *i,
                j: *j,
                eps: *si,
                tau: *sj,
            },
            _ => unreachable!(),
        };
        expect.push(Expect {
            target: k,
            rel: want,
            kind: "production",
        });
        n_prod += 1;
    }
    println!(
        "  targets: {} positive + 16 degenerate controls, {} random, {} production (collect_relations: {} targets scanned, {} relations)",
        n_controls, n_random, n_prod, stats.targets, stats.relations
    );

    // ── CPU reference ────────────────────────────────────────────
    let t1 = Instant::now();
    let reference = cpu_reference64(&curve, &fb, &targets);
    let cpu_secs = t1.elapsed().as_secs_f64();
    println!(
        "  cpu prime_fast reference: {} relations over {} targets in {:.3}s ({:.3e} (target,i) pairs/s)",
        reference.len(),
        targets.len(),
        cpu_secs,
        (targets.len() * fb.len()) as f64 / cpu_secs
    );
    let mut failures = 0usize;
    for r in &reference {
        if !verify_relation64(&curve, &fb, targets[r.target as usize], r) {
            println!("  FAIL reference relation does not re-add: {r:?}");
            failures += 1;
        }
    }

    // ── oracle, every configuration ──────────────────────────────
    let (fbx, fby) = factor_base_canonical(&curve, &fb);
    let (tx, ty) = targets_canonical(&curve, &targets);
    let mut rows = Vec::new();
    for lookup in [Lookup::Sorted, Lookup::Hash] {
        for chunk in [16u32, 32, 64] {
            let opts = Options {
                backend,
                lookup,
                chunk,
            };
            let t2 = Instant::now();
            let rep = match decompose64(curve.f.p, &fbx, &fby, &tx, &ty, opts) {
                Ok(r) => r,
                Err(e) => {
                    println!("  FAIL {lookup:?} chunk={chunk}: {e}");
                    failures += 1;
                    continue;
                }
            };
            let wall = t2.elapsed().as_secs_f64();
            let got: BTreeSet<Relation> = rep.relations.iter().copied().collect();
            let mut bad_readd = 0usize;
            for r in &rep.relations {
                if !verify_relation64(&curve, &fb, targets[r.target as usize], r) {
                    bad_readd += 1;
                }
            }
            let dup = rep.relations.len() - got.len();
            let missing: Vec<&Relation> = reference.difference(&got).collect();
            let spurious: Vec<&Relation> = got.difference(&reference).collect();
            let mut missing_ctl = 0usize;
            for e in &expect {
                if !got.contains(&e.rel) {
                    missing_ctl += 1;
                    if missing_ctl <= 5 {
                        println!(
                            "  FAIL {lookup:?} chunk={chunk}: {} control target {} missing {:?}",
                            e.kind, e.target, e.rel
                        );
                    }
                }
            }
            let ok = bad_readd == 0
                && dup == 0
                && missing.is_empty()
                && spurious.is_empty()
                && missing_ctl == 0;
            if !ok {
                failures += 1;
                for r in missing.iter().take(5) {
                    println!("  FAIL missing vs cpu: {r:?}");
                }
                for r in spurious.iter().take(5) {
                    println!("  FAIL spurious vs cpu: {r:?}");
                }
            }
            let pairs = (targets.len() * fb.len()) as f64;
            println!(
                "  {:<7} {:<6} chunk={:<2} : {}  {} relations, re-add fail {}, dup {}, missing {}, spurious {}, controls missing {}  [{:.3}s kernel, {:.3}s wall, {:.3e} pairs/s]",
                format!("{:?}", rep.backend_used).to_lowercase(),
                format!("{:?}", rep.lookup_used).to_lowercase(),
                rep.chunk_used,
                if ok { "PASS" } else { "FAIL" },
                rep.relations.len(),
                bad_readd,
                dup,
                missing.len(),
                spurious.len(),
                missing_ctl,
                rep.seconds,
                wall,
                pairs / rep.seconds.max(1e-9)
            );
            rows.push(serde_json::json!({
                "backend": format!("{:?}", rep.backend_used).to_lowercase(),
                "lookup": format!("{:?}", rep.lookup_used).to_lowercase(),
                "chunk": rep.chunk_used,
                "pass": ok,
                "relations": rep.relations.len(),
                "readd_failures": bad_readd,
                "duplicates": dup,
                "missing_vs_cpu": missing.len(),
                "spurious_vs_cpu": spurious.len(),
                "controls_missing": missing_ctl,
                "kernel_seconds": rep.seconds,
                "wall_seconds": wall,
                "pairs_per_second": pairs / rep.seconds.max(1e-9),
            }));
        }
    }

    let verdict = if failures == 0 { "PASS" } else { "FAIL" };
    println!(
        "RESULT {verdict}: {} targets ({} controls, {} random, {} production), {} reference relations, {} configurations, {} failures",
        targets.len(),
        n_ctl,
        n_random,
        n_prod,
        reference.len(),
        rows.len(),
        failures
    );

    if let Ok(dir) = std::env::var("BENCH_OUT") {
        std::fs::create_dir_all(&dir).expect("BENCH_OUT dir");
        let doc = serde_json::json!({
            "example": "gpu_rootfind_validate",
            "date": "2026-10-12",
            "backend_library": backend_info(),
            "cuda_available": cuda_available(),
            "curve": { "name": curve.name, "p": curve.f.p.to_string(), "bits": bits, "n": curve.n.to_string(), "a": "-3" },
            "factor_base": fb.len(),
            "targets": { "total": targets.len(), "positive_controls": n_controls, "degenerate_controls": 16, "random": n_random, "production": n_prod },
            "collect_relations": { "targets_scanned": stats.targets, "relations": stats.relations },
            "cpu_reference": { "relations": reference.len(), "seconds": cpu_secs },
            "configurations": rows,
            "failures": failures,
            "verdict": verdict,
            "claim_boundary": "decomposition subroutine on public synthetic known-answer instances; emulation backend; no key-recovery, sub-rho or deployed-curve claim",
        });
        let path = format!("{dir}/validate_{bits}bit_fb{}_seed{seed}.json", fb.len());
        std::fs::write(&path, serde_json::to_string_pretty(&doc).unwrap()).expect("write json");
        println!("wrote {path}");
    }
    if failures != 0 {
        std::process::exit(1);
    }
}
