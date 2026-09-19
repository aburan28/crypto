//! Speed harness for the post-quantum schemes.
//!
//! One table, one unit: **cycles per operation**, median over repetitions,
//! with a correctness column.  The unit is cycles rather than wall time so
//! that the numbers can be compared against the published figures for the
//! reference implementations, which are all quoted in cycles.
//!
//! The reference column is the boundary in the sense of AGENTS.md: the cost of
//! the best implementation that already does the same job, in the same unit.
//! For ML-KEM and ML-DSA those are the pq-crystals C reference and AVX2
//! implementations; the numbers are recorded in `docs/pqc-speed.md` with their
//! source, not guessed here.
//!
//! Run with:
//!
//! ```sh
//! cargo bench --bench pqc_speed
//! ```
//!
//! or, for a single scheme, `cargo bench --bench pqc_speed -- ml-kem`.

use std::time::Instant;

// ── cycle counter ────────────────────────────────────────────────────────────

/// Read the cycle counter.  On x86_64 this is `rdtsc`, which counts at a fixed
/// reference frequency rather than the core clock; on everything else it is a
/// nanosecond timer scaled by the measured reference frequency.  Either way the
/// unit is consistent within a run, which is what the table needs.
#[cfg(target_arch = "x86_64")]
#[inline]
fn cycles() -> u64 {
    // SAFETY: _rdtsc is available on every x86_64 CPU (Pentium and later) and
    // has no memory operands or side effects.
    unsafe { core::arch::x86_64::_rdtsc() }
}

#[cfg(not(target_arch = "x86_64"))]
#[inline]
fn cycles() -> u64 {
    use std::time::SystemTime;
    SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap()
        .as_nanos() as u64
}

/// Cycles per second of whatever `cycles()` counts, measured once.
fn cycle_hz() -> f64 {
    let t0 = Instant::now();
    let c0 = cycles();
    while t0.elapsed().as_millis() < 200 {
        std::hint::spin_loop();
    }
    let c1 = cycles();
    let dt = t0.elapsed().as_secs_f64();
    (c1 - c0) as f64 / dt
}

// ── measurement ──────────────────────────────────────────────────────────────

/// Median cycles for one call of `f`, over `reps` repetitions after `warmup`.
///
/// The median, not the mean: a scheduler preemption in the middle of a run adds
/// a six-figure outlier that would move a mean by more than any optimisation.
fn measure<T>(warmup: usize, reps: usize, mut f: impl FnMut() -> T) -> f64 {
    for _ in 0..warmup {
        std::hint::black_box(f());
    }
    let mut samples = Vec::with_capacity(reps);
    for _ in 0..reps {
        let c0 = cycles();
        std::hint::black_box(f());
        let c1 = cycles();
        samples.push(c1.wrapping_sub(c0));
    }
    samples.sort_unstable();
    samples[samples.len() / 2] as f64
}

/// Median cycles for *one* operation, where `f` performs `inner` of them.
///
/// `rdtsc` costs on the order of the thing being measured once an operation is
/// down at tens of cycles, so the field primitives are timed in batches and
/// divided.  The batches are dependency chains (`x = x * b`), so these are
/// latency numbers, not throughput: the pessimistic end of the range.
fn measure_each<T>(warmup: usize, reps: usize, inner: usize, f: impl FnMut() -> T) -> f64 {
    measure(warmup, reps, f) / inner as f64
}

struct Row {
    scheme: &'static str,
    op: &'static str,
    cycles: f64,
    ok: bool,
}

fn print_table(rows: &[Row], hz: f64) {
    println!();
    println!(
        "| scheme | operation | kilocycles | ops/s | correct |\n\
         |--------|-----------|-----------:|------:|:-------:|"
    );
    for r in rows {
        println!(
            "| {} | {} | {:.3} | {:.0} | {} |",
            r.scheme,
            r.op,
            r.cycles / 1000.0,
            hz / r.cycles,
            if r.ok { "yes" } else { "NO" }
        );
    }
    println!();
}

// ── ML-KEM ───────────────────────────────────────────────────────────────────

fn bench_ml_kem_fast(rows: &mut Vec<Row>) {
    use crypto_lib::pqc::fast::ml_kem as fast;
    use crypto_lib::pqc::ml_kem::{ML_KEM_1024, ML_KEM_512, ML_KEM_768};

    for (name, params) in [
        ("ML-KEM-512 fast", ML_KEM_512),
        ("ML-KEM-768 fast", ML_KEM_768),
        ("ML-KEM-1024 fast", ML_KEM_1024),
    ] {
        let p = &params;
        let (ek, dk) = fast::ml_kem_keygen(p);
        let (ct, k_enc) = fast::ml_kem_encaps(p, &ek).expect("encaps");
        let k_dec = fast::ml_kem_decaps(p, &dk, &ct).expect("decaps");
        // The correctness column also holds the fast output against the
        // reference decapsulation, so a wrong-but-self-consistent
        // implementation cannot show up here as a speedup.
        let agrees = crypto_lib::pqc::ml_kem::ml_kem_decaps(p, &dk, &ct) == Some(k_enc);
        let ok = k_enc == k_dec && agrees;

        rows.push(Row {
            scheme: name,
            op: "keygen",
            cycles: measure(3, 25, || fast::ml_kem_keygen(p)),
            ok,
        });
        rows.push(Row {
            scheme: name,
            op: "encaps",
            cycles: measure(3, 25, || fast::ml_kem_encaps(p, &ek).unwrap()),
            ok,
        });
        rows.push(Row {
            scheme: name,
            op: "decaps",
            cycles: measure(3, 25, || fast::ml_kem_decaps(p, &dk, &ct).unwrap()),
            ok,
        });
    }
}

fn bench_ml_kem(rows: &mut Vec<Row>) {
    use crypto_lib::pqc::ml_kem::*;

    for (name, params) in [
        ("ML-KEM-512", ML_KEM_512),
        ("ML-KEM-768", ML_KEM_768),
        ("ML-KEM-1024", ML_KEM_1024),
    ] {
        let p = &params;
        let (ek, dk) = ml_kem_keygen(p);

        // One round trip up front: the correctness column is the shared secret
        // agreeing, checked once per parameter set rather than per sample.
        let (ct, k_enc) = ml_kem_encaps(p, &ek).expect("encaps");
        let k_dec = ml_kem_decaps(p, &dk, &ct).expect("decaps");
        let ok = k_enc == k_dec;

        rows.push(Row {
            scheme: name,
            op: "keygen",
            cycles: measure(3, 25, || ml_kem_keygen(p)),
            ok,
        });
        rows.push(Row {
            scheme: name,
            op: "encaps",
            cycles: measure(3, 25, || ml_kem_encaps(p, &ek).unwrap()),
            ok,
        });
        rows.push(Row {
            scheme: name,
            op: "decaps",
            cycles: measure(3, 25, || ml_kem_decaps(p, &dk, &ct).unwrap()),
            ok,
        });
    }
}

// ── ML-DSA ───────────────────────────────────────────────────────────────────

fn bench_ml_dsa(rows: &mut Vec<Row>) {
    use crypto_lib::pqc::ml_dsa::*;

    let seed = [7u8; 32];
    let rnd = [0u8; 32];
    let msg = b"the quick brown fox jumps over the lazy dog";

    let (pk, sk) = ml_dsa_65_keygen(&seed);
    let sig = ml_dsa_65_sign(&sk, msg, &rnd);
    let ok = ml_dsa_65_verify(&pk, msg, &sig);

    rows.push(Row {
        scheme: "ML-DSA-65",
        op: "keygen",
        cycles: measure(2, 11, || ml_dsa_65_keygen(&seed)),
        ok,
    });
    rows.push(Row {
        scheme: "ML-DSA-65",
        op: "sign",
        cycles: measure(2, 11, || ml_dsa_65_sign(&sk, msg, &rnd)),
        ok,
    });
    rows.push(Row {
        scheme: "ML-DSA-65",
        op: "verify",
        cycles: measure(2, 11, || ml_dsa_65_verify(&pk, msg, &sig)),
        ok,
    });
}

fn bench_ml_dsa_fast(rows: &mut Vec<Row>) {
    use crypto_lib::pqc::fast::ml_dsa as fast;
    use crypto_lib::pqc::ml_dsa as reference;

    let seed = [7u8; 32];
    let rnd = [0u8; 32];
    let msg = b"the quick brown fox jumps over the lazy dog";

    let (pk, sk) = fast::ml_dsa_65_keygen(&seed);
    let sig = fast::ml_dsa_65_sign(&sk, msg, &rnd);

    // The correctness column is the round trip *and* byte equality with the
    // reference, so a wrong-but-self-consistent implementation cannot appear
    // here as a speedup.
    let (rpk, rsk) = reference::ml_dsa_65_keygen(&seed);
    let rsig = reference::ml_dsa_65_sign(&rsk, msg, &rnd);
    let ok = fast::ml_dsa_65_verify(&pk, msg, &sig)
        && pk.0 == rpk.0
        && sk.0 == rsk.0
        && sig == rsig
        && reference::ml_dsa_65_verify(&rpk, msg, &sig)
        && fast::ml_dsa_65_verify(&pk, msg, &rsig);

    rows.push(Row {
        scheme: "ML-DSA-65 fast",
        op: "keygen",
        cycles: measure(2, 11, || fast::ml_dsa_65_keygen(&seed)),
        ok,
    });
    rows.push(Row {
        scheme: "ML-DSA-65 fast",
        op: "sign",
        cycles: measure(2, 11, || fast::ml_dsa_65_sign(&sk, msg, &rnd)),
        ok,
    });
    rows.push(Row {
        scheme: "ML-DSA-65 fast",
        op: "verify",
        cycles: measure(2, 11, || fast::ml_dsa_65_verify(&pk, msg, &sig)),
        ok,
    });
}

// ── SQIsign ──────────────────────────────────────────────────────────────────

fn bench_sqisign(rows: &mut Vec<Row>) {
    use crypto_lib::pqc::sqisign::*;

    let msg = b"the quick brown fox jumps over the lazy dog";
    let (pk, sk) = sqisign_keygen();
    let sig = sqisign_sign(&pk, &sk, msg);
    let ok = sqisign_verify(&pk, msg, &sig);

    rows.push(Row {
        scheme: "SQIsign (toy p=431)",
        op: "keygen",
        cycles: measure(2, 11, sqisign_keygen),
        ok,
    });
    rows.push(Row {
        scheme: "SQIsign (toy p=431)",
        op: "sign",
        cycles: measure(2, 11, || sqisign_sign(&pk, &sk, msg)),
        ok,
    });
    rows.push(Row {
        scheme: "SQIsign (toy p=431)",
        op: "verify",
        cycles: measure(2, 11, || sqisign_verify(&pk, msg, &sig)),
        ok,
    });
}

// ── primitives ───────────────────────────────────────────────────────────────

/// The hash layer, separately: in the reference implementations of both lattice
/// schemes Keccak is the single largest line item, so its cost belongs in the
/// table on its own before anything above it is attributed.
fn bench_primitives(rows: &mut Vec<Row>) {
    use crypto_lib::hash::sha3::{sha3_256, sha3_512, shake128, shake256};

    let block = [0u8; 32];
    rows.push(Row {
        scheme: "primitive",
        op: "sha3-256(32B)",
        cycles: measure(50, 201, || sha3_256(&block)),
        ok: true,
    });
    rows.push(Row {
        scheme: "primitive",
        op: "sha3-512(32B)",
        cycles: measure(50, 201, || sha3_512(&block)),
        ok: true,
    });
    rows.push(Row {
        scheme: "primitive",
        op: "shake128(34B->504B)",
        cycles: measure(50, 201, || shake128(&block, 504)),
        ok: true,
    });
    rows.push(Row {
        scheme: "primitive",
        op: "shake256(32B->128B)",
        cycles: measure(50, 201, || shake256(&block, 128)),
        ok: true,
    });
}

fn main() {
    let filter: Vec<String> = std::env::args().skip(1).filter(|a| !a.starts_with('-')).collect();
    let want = |name: &str| filter.is_empty() || filter.iter().any(|f| name.contains(f.as_str()));

    let hz = cycle_hz();
    println!("reference clock: {:.3} GHz", hz / 1e9);

    let mut rows = Vec::new();
    if want("primitives") {
        bench_primitives(&mut rows);
    }
    if want("ml-kem") {
        bench_ml_kem(&mut rows);
        bench_ml_kem_fast(&mut rows);
    }
    if want("ml-dsa") {
        bench_ml_dsa(&mut rows);
        bench_ml_dsa_fast(&mut rows);
    }
    if want("sqisign") {
        bench_sqisign(&mut rows);
    }
    print_table(&rows, hz);

    if rows.iter().any(|r| !r.ok) {
        eprintln!("error: a scheme failed its round trip; the numbers above are meaningless");
        std::process::exit(1);
    }
}
