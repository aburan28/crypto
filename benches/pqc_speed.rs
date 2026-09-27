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

// ── SQIsign field / isogeny core ─────────────────────────────────────────────

/// The layer SQIsign *verification* actually spends its time in, at the real
/// NIST level-1 parameter size: `p = 3*2^324 - 1`, `F_p^2`, x-only Montgomery
/// curves and 2-isogeny chains.  See `src/pqc/fast/isogeny.rs`; this is the
/// arithmetic core only, not the scheme.
fn bench_isogeny(rows: &mut Vec<Row>) {
    use crypto_lib::pqc::fast::isogeny::*;

    const S: &str = "isogeny p=3*2^324-1";

    // ── field ────────────────────────────────────────────────────────────
    let a = Fp::from_u64(0x0123_4567_89ab_cdef).mul(&Fp::from_u64(0xfedc_ba98_7654_3210));
    let b = Fp::from_u64(0x9e37_79b9_7f4a_7c15);
    let a2 = Fp2::new(a, b);
    let b2 = Fp2::new(b, a);

    let field_ok = a.mul(&a.inv()) == Fp::ONE && a2.mul(&a2.inv()) == Fp2::ONE;

    rows.push(Row {
        scheme: S,
        op: "F_p mul (6 limbs)",
        cycles: measure_each(200, 401, 128, || {
            let mut x = std::hint::black_box(a);
            for _ in 0..128 {
                x = x.mul(&b);
            }
            x
        }),
        ok: field_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "F_p sqr (= mul, see docs)",
        cycles: measure_each(200, 401, 128, || {
            let mut x = std::hint::black_box(a);
            for _ in 0..128 {
                x = x.sqr();
            }
            x
        }),
        ok: field_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "F_p inv (addition chain)",
        cycles: measure_each(5, 51, 4, || {
            let mut x = std::hint::black_box(a);
            for _ in 0..4 {
                x = x.inv();
            }
            x
        }),
        ok: field_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "F_p^2 mul (karatsuba)",
        cycles: measure_each(200, 401, 128, || {
            let mut x = std::hint::black_box(a2);
            for _ in 0..128 {
                x = x.mul(&b2);
            }
            x
        }),
        ok: field_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "F_p^2 sqr",
        cycles: measure_each(200, 401, 128, || {
            let mut x = std::hint::black_box(a2);
            for _ in 0..128 {
                x = x.sqr();
            }
            x
        }),
        ok: field_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "F_p^2 inv",
        cycles: measure_each(5, 51, 4, || {
            let mut x = std::hint::black_box(a2);
            for _ in 0..4 {
                x = x.inv();
            }
            x
        }),
        ok: field_ok,
    });

    // ── curve ────────────────────────────────────────────────────────────
    // E_6: y^2 = x^3 + 6x^2 + x, supersingular over F_p^2 with
    // #E = (p+1)^2 = (3*2^324)^2.  Find a point of E by trying small x, which
    // is deterministic and needs no RNG in the harness.
    let curve = Curve::from_a(Fp2::from_u64(6));
    let a24 = curve.normalised_a24();

    // Enumerate x = c + d*i for small c, d, which is deterministic and needs
    // no RNG in the harness.  d must be nonzero: E_6 has its full 2-torsion
    // rational over F_p (2 is a QR mod p since p = 7 mod 8), so E(F_p) has
    // 2-Sylow Z/2 x Z/2^323 and no x in F_p carries a point of order 2^324.
    let candidates = (1u64..40).flat_map(|d| (0u64..40).map(move |c| (c, d)));

    let mut base = PointX::INFINITY;
    let mut kernel = PointX::INFINITY;
    for (c, d) in candidates {
        let x = Fp2::new(Fp::from_u64(c), Fp::from_u64(d));
        let p = PointX::from_affine(x);
        if !curve.is_on_curve(&p) {
            continue;
        }
        if base.is_infinity() {
            base = p;
        }
        // A point of exact order 2^324 whose bottom 2-torsion point is not
        // (0,0), the kernel the x-only isogeny formulas exclude.
        let q = ladder(&[3], 8, &p, &a24);
        if q.is_infinity() {
            continue;
        }
        let bottom = xdbl_e(&q, TWO_TORSION_POWER - 1, &a24);
        if bottom.is_infinity() || bottom.x.is_zero() {
            continue;
        }
        kernel = q;
        break;
    }
    assert!(!base.is_infinity(), "no base point found");
    assert!(!kernel.is_infinity(), "no usable full-order 2-torsion point found");
    let curve_ok = xdbl_e(&kernel, TWO_TORSION_POWER, &a24).is_infinity()
        && ladder(&P_PLUS_1, P_PLUS_1_BITS, &base, &a24).is_infinity();

    rows.push(Row {
        scheme: S,
        op: "xDBL (A24plus:C24)",
        cycles: measure_each(100, 301, 64, || {
            let mut p = std::hint::black_box(base);
            for _ in 0..64 {
                p = xdbl(&p, &curve);
            }
            p
        }),
        ok: curve_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "xDBL (a24 normalised)",
        cycles: measure_each(100, 301, 64, || {
            let mut p = std::hint::black_box(base);
            for _ in 0..64 {
                p = xdbl_a24(&p, &a24);
            }
            p
        }),
        ok: curve_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "xADD",
        cycles: measure_each(100, 301, 64, || {
            let d = xdbl(&base, &curve);
            let mut p = std::hint::black_box(d);
            for _ in 0..64 {
                p = xadd(&p, &base, &base);
            }
            p
        }),
        ok: curve_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "ladder [k]P, k ~ 2^326",
        cycles: measure(3, 31, || ladder(&P_PLUS_1, P_PLUS_1_BITS, &base, &a24)),
        ok: curve_ok,
    });

    // ── isogenies ────────────────────────────────────────────────────────
    let k2 = xdbl_e(&kernel, TWO_TORSION_POWER - 1, &a24); // order 2
    let k4 = xdbl_e(&kernel, TWO_TORSION_POWER - 2, &a24); // order 4
    let (_, kps2) = isog2_codomain(&k2).expect("generic kernel");
    let (_, kps4) = isog4_codomain(&k4);
    let isog_ok = isog2_eval(&k2, &kps2).is_infinity() && isog4_eval(&k4, &kps4).is_infinity();

    rows.push(Row {
        scheme: S,
        op: "2-isog step (codomain+eval)",
        cycles: measure_each(100, 301, 64, || {
            let k = std::hint::black_box(k2);
            let mut q = std::hint::black_box(base);
            for _ in 0..64 {
                let (_, kps) = isog2_codomain(&k).unwrap();
                q = isog2_eval(&q, &kps);
            }
            q
        }),
        ok: isog_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "2-isog eval only",
        cycles: measure_each(100, 301, 64, || {
            let mut q = std::hint::black_box(base);
            for _ in 0..64 {
                q = isog2_eval(&q, &kps2);
            }
            q
        }),
        ok: isog_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "4-isog step (codomain+eval)",
        cycles: measure_each(100, 301, 64, || {
            let k = std::hint::black_box(k4);
            let mut q = std::hint::black_box(base);
            for _ in 0..64 {
                let (_, kps) = isog4_codomain(&k);
                q = isog4_eval(&q, &kps);
            }
            q
        }),
        ok: isog_ok,
    });

    // ── full chains ──────────────────────────────────────────────────────
    let n = TWO_TORSION_POWER / 2; // 162 four-isogeny steps = degree 2^324
    let strategy = optimal_strategy(n, COST_DBL, COST_EVAL);

    let mut probe = [kernel, base];
    let chain_ok = two_isogeny_chain_with_strategy(
        &curve,
        &kernel,
        TWO_TORSION_POWER,
        &mut probe,
        &strategy,
    )
    .map(|img| probe[0].is_infinity() && img.is_on_curve(&probe[1]))
    .unwrap_or(false);

    rows.push(Row {
        scheme: S,
        op: "2^324 chain, optimal strategy",
        cycles: measure(2, 11, || {
            let mut push = [base];
            two_isogeny_chain_with_strategy(
                &curve,
                &kernel,
                TWO_TORSION_POWER,
                &mut push,
                &strategy,
            )
        }),
        ok: chain_ok,
    });
    rows.push(Row {
        scheme: S,
        op: "2^324 chain, naive walker",
        cycles: measure(1, 5, || {
            let mut push = [base];
            chain_4_naive(&curve, &kernel, n, &mut push)
        }),
        ok: chain_ok,
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
    if want("isogeny") {
        bench_isogeny(&mut rows);
    }
    print_table(&rows, hz);

    if rows.iter().any(|r| !r.ok) {
        eprintln!("error: a scheme failed its round trip; the numbers above are meaningless");
        std::process::exit(1);
    }
}
