//! `crypto pqc-fast …` — the command-line surface for
//! [`crypto_lib::pqc::fast`].
//!
//! The speed-oriented ML-KEM, ML-DSA, Keccak and isogeny code had no entry
//! point outside the benchmark harness and its own tests. It does now: every
//! module here is reachable, self-testable and benchmarkable from the CLI.
//!
//! The library's warning applies with full force to everything this command
//! runs: **none of it is constant-time**, and the attacks under
//! `crypto mlwe kem-pco` and `crypto mlwe dsa-leak` are exactly what that costs.

use clap::Subcommand;
use crypto_lib::pqc::fast::isogeny as iso;
use crypto_lib::pqc::fast::keccak;
use crypto_lib::pqc::fast::ml_dsa as fast_dsa;
use crypto_lib::pqc::fast::ml_kem as fast_kem;
use crypto_lib::pqc::ml_kem::{MlKemParams, ML_KEM_1024, ML_KEM_512, ML_KEM_768};
use crypto_lib::utils::encoding::to_hex;
use std::time::Instant;

#[derive(Subcommand)]
pub enum PqcFastOp {
    /// The speed-oriented ML-KEM: keygen, encapsulate, decapsulate.
    MlKem {
        /// `512`, `768` or `1024`.
        #[arg(long, default_value = "768")]
        param: String,
        /// `demo`, `selftest` or `bench`.
        #[arg(long, default_value = "demo")]
        op: String,
        #[arg(long, default_value_t = 200)]
        iters: u32,
    },
    /// The speed-oriented ML-DSA-65: keygen, sign, verify.
    MlDsa {
        /// `demo`, `selftest` or `bench`.
        #[arg(long, default_value = "demo")]
        op: String,
        #[arg(long, default_value = "post-quantum signatures")]
        message: String,
        #[arg(long, default_value_t = 50)]
        iters: u32,
    },
    /// The unrolled Keccak-f[1600] and its sponges.
    Keccak {
        /// `demo`, `selftest` or `bench`.
        #[arg(long, default_value = "demo")]
        op: String,
        #[arg(long, default_value = "")]
        message: String,
        #[arg(long, default_value_t = 2000)]
        iters: u32,
    },
    /// The SQIsign-scale isogeny arithmetic: F_p² over a 326-bit prime,
    /// x-only Montgomery curves, and 2-power isogeny chains.
    Isogeny {
        /// `info`, `selftest`, `chain`, `strategy` or `bench`.
        #[arg(long, default_value = "info")]
        op: String,
        /// Chain length for `chain`; the SQIsign level-1 value is 324.
        #[arg(long, default_value_t = 324)]
        e: usize,
        #[arg(long, default_value_t = 200)]
        iters: u32,
    },
    /// Run every module's self-test.
    Selftest,
    /// Run every module's benchmark.
    Bench,
}

fn params(name: &str) -> &'static MlKemParams {
    let key: String = name.chars().filter(|c| c.is_ascii_digit()).collect();
    match key.as_str() {
        "512" => &ML_KEM_512,
        "1024" => &ML_KEM_1024,
        _ => &ML_KEM_768,
    }
}

/// Time `iters` calls and return microseconds per call.
fn timed<F: FnMut()>(iters: u32, mut f: F) -> f64 {
    // One warm-up call, so the first allocation and any lazy table setup is
    // not charged to the measurement.
    f();
    let t = Instant::now();
    for _ in 0..iters {
        f();
    }
    t.elapsed().as_secs_f64() * 1e6 / iters as f64
}

pub fn run(op: PqcFastOp) {
    match op {
        PqcFastOp::MlKem { param, op, iters } => ml_kem(params(&param), &op, iters),
        PqcFastOp::MlDsa { op, message, iters } => ml_dsa(&op, message.as_bytes(), iters),
        PqcFastOp::Keccak { op, message, iters } => keccak_cmd(&op, message.as_bytes(), iters),
        PqcFastOp::Isogeny { op, e, iters } => isogeny(&op, e, iters),
        PqcFastOp::Selftest => {
            ml_kem(&ML_KEM_512, "selftest", 0);
            ml_kem(&ML_KEM_768, "selftest", 0);
            ml_kem(&ML_KEM_1024, "selftest", 0);
            ml_dsa("selftest", b"selftest", 0);
            keccak_cmd("selftest", b"", 0);
            isogeny("selftest", 324, 0);
        }
        PqcFastOp::Bench => {
            ml_kem(&ML_KEM_768, "bench", 200);
            ml_dsa("bench", b"benchmark", 50);
            keccak_cmd("bench", b"", 2000);
            isogeny("bench", 324, 200);
        }
    }
}

fn ml_kem(p: &MlKemParams, op: &str, iters: u32) {
    match op {
        "demo" => {
            let (ek, dk) = fast_kem::ml_kem_keygen(p);
            println!("{} (pqc::fast)", p.name);
            println!("  ek: {} bytes, dk: {} bytes", ek.0.len(), dk.0.len());
            let (ct, ss) = fast_kem::ml_kem_encaps(p, &ek).expect("fresh key is valid");
            println!("  ct: {} bytes", ct.len());
            println!("  shared secret (encaps): {}", to_hex(&ss));
            let ss2 = fast_kem::ml_kem_decaps(p, &dk, &ct).expect("well-formed ciphertext");
            println!("  shared secret (decaps): {}", to_hex(&ss2));
            println!("  match: {}", ss == ss2);
        }
        "selftest" => {
            let mut ok = true;
            // Round trip.
            let (ek, dk) = fast_kem::ml_kem_keygen_internal(p, &[7u8; 32], &[9u8; 32]);
            let (ct, ss) = fast_kem::ml_kem_encaps_internal(p, &ek, &[3u8; 32]).unwrap();
            ok &= fast_kem::ml_kem_decaps(p, &dk, &ct) == Some(ss);
            // Agreement with the reference implementation, byte for byte.
            let (rek, rdk) =
                crypto_lib::pqc::ml_kem::ml_kem_keygen_internal(p, &[7u8; 32], &[9u8; 32]);
            ok &= rek.0 == ek.0 && rdk.0 == dk.0;
            let (rct, rss) =
                crypto_lib::pqc::ml_kem::ml_kem_encaps_internal(p, &rek, &[3u8; 32]).unwrap();
            ok &= rct == ct && rss == ss;
            // Implicit rejection: a tampered ciphertext must not error, and
            // must not yield the honest secret.
            let mut bad = ct.clone();
            bad[0] ^= 1;
            let rej = fast_kem::ml_kem_decaps(p, &dk, &bad);
            ok &= rej.is_some() && rej != Some(ss);
            // Key check accepts the real key and rejects a mangled one.
            ok &= fast_kem::ml_kem_check_ek(p, &ek.0);
            let mut bad_ek = ek.0.clone();
            bad_ek[0] = 0xff;
            bad_ek[1] = 0xff;
            ok &= !fast_kem::ml_kem_check_ek(p, &bad_ek);
            println!("{} selftest: {}", p.name, if ok { "PASS" } else { "FAIL" });
            if !ok {
                std::process::exit(1);
            }
        }
        "bench" => {
            let (ek, dk) = fast_kem::ml_kem_keygen(p);
            let (ct, _) = fast_kem::ml_kem_encaps(p, &ek).unwrap();
            let kg = timed(iters, || {
                std::hint::black_box(fast_kem::ml_kem_keygen_internal(p, &[1u8; 32], &[2u8; 32]));
            });
            let en = timed(iters, || {
                std::hint::black_box(fast_kem::ml_kem_encaps_internal(p, &ek, &[3u8; 32]));
            });
            let de = timed(iters, || {
                std::hint::black_box(fast_kem::ml_kem_decaps(p, &dk, &ct));
            });
            println!("{} (pqc::fast), {iters} iterations", p.name);
            println!("  keygen:  {kg:>9.2} us");
            println!("  encaps:  {en:>9.2} us");
            println!("  decaps:  {de:>9.2} us");
        }
        other => {
            eprintln!("unknown op `{other}`; try demo, selftest or bench");
            std::process::exit(2);
        }
    }
}

fn ml_dsa(op: &str, msg: &[u8], iters: u32) {
    match op {
        "demo" => {
            let (pk, sk) = fast_dsa::ml_dsa_65_keygen(&[5u8; 32]);
            let sig = fast_dsa::ml_dsa_65_sign(&sk, msg, &[0u8; 32]);
            println!("ML-DSA-65 (pqc::fast)");
            println!("  pk: {} bytes, sk: {} bytes", pk.0.len(), sk.0.len());
            println!("  message: {} bytes", msg.len());
            println!("  signature: {} bytes", sig.len());
            println!("  verify: {}", fast_dsa::ml_dsa_65_verify(&pk, msg, &sig));
        }
        "selftest" => {
            let mut ok = true;
            let (pk, sk) = fast_dsa::ml_dsa_65_keygen(&[5u8; 32]);
            let sig = fast_dsa::ml_dsa_65_sign(&sk, msg, &[0u8; 32]);
            ok &= fast_dsa::ml_dsa_65_verify(&pk, msg, &sig);
            // A tampered signature must be rejected.
            let mut bad = sig.clone();
            bad[0] ^= 1;
            ok &= !fast_dsa::ml_dsa_65_verify(&pk, msg, &bad);
            // A different message must be rejected.
            ok &= !fast_dsa::ml_dsa_65_verify(&pk, b"a different message", &sig);
            // Byte-for-byte agreement with the reference implementation.
            let (rpk, rsk) = crypto_lib::pqc::ml_dsa::ml_dsa_65_keygen(&[5u8; 32]);
            ok &= rpk.0 == pk.0 && rsk.0 == sk.0;
            let rsig = crypto_lib::pqc::ml_dsa::ml_dsa_65_sign(&rsk, msg, &[0u8; 32]);
            ok &= rsig == sig;
            // Cross-verification both ways.
            ok &= crypto_lib::pqc::ml_dsa::ml_dsa_65_verify(&rpk, msg, &sig);
            ok &= fast_dsa::ml_dsa_65_verify(&pk, msg, &rsig);
            println!("ML-DSA-65 selftest: {}", if ok { "PASS" } else { "FAIL" });
            if !ok {
                std::process::exit(1);
            }
        }
        "bench" => {
            let (pk, sk) = fast_dsa::ml_dsa_65_keygen(&[5u8; 32]);
            let sig = fast_dsa::ml_dsa_65_sign(&sk, msg, &[0u8; 32]);
            let kg = timed(iters, || {
                std::hint::black_box(fast_dsa::ml_dsa_65_keygen(&[5u8; 32]));
            });
            let sg = timed(iters, || {
                std::hint::black_box(fast_dsa::ml_dsa_65_sign(&sk, msg, &[0u8; 32]));
            });
            let vf = timed(iters, || {
                std::hint::black_box(fast_dsa::ml_dsa_65_verify(&pk, msg, &sig));
            });
            println!("ML-DSA-65 (pqc::fast), {iters} iterations");
            println!("  keygen:  {kg:>9.2} us");
            println!("  sign:    {sg:>9.2} us");
            println!("  verify:  {vf:>9.2} us");
        }
        other => {
            eprintln!("unknown op `{other}`; try demo, selftest or bench");
            std::process::exit(2);
        }
    }
}

fn keccak_cmd(op: &str, msg: &[u8], iters: u32) {
    match op {
        "demo" => {
            let mut out = [0u8; 32];
            keccak::shake128_into(&mut out, msg);
            println!("Keccak (pqc::fast)");
            println!("  SHAKE128(msg, 32): {}", to_hex(&out));
            keccak::shake256_into(&mut out, msg);
            println!("  SHAKE256(msg, 32): {}", to_hex(&out));
            println!("  SHA3-256(msg):     {}", to_hex(&keccak::sha3_256(msg)));
        }
        "selftest" => {
            let mut ok = true;
            // Against the library's own reference SHA-3, on several lengths —
            // including ones that straddle a rate boundary.
            for len in [0usize, 1, 135, 136, 137, 271, 1000] {
                let m: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
                ok &= keccak::sha3_256(&m) == crypto_lib::hash::sha3::sha3_256(&m);
                let mut a = vec![0u8; 97];
                keccak::shake128_into(&mut a, &m);
                ok &= a == crypto_lib::hash::sha3::shake128(&m, 97);
                let mut b = vec![0u8; 97];
                keccak::shake256_into(&mut b, &m);
                ok &= b == crypto_lib::hash::sha3::shake256(&m, 97);
            }
            // The two-part helpers must match the concatenated one-part call.
            let (x, y) = (b"left".as_slice(), b"right".as_slice());
            let mut joined = x.to_vec();
            joined.extend_from_slice(y);
            let mut two = [0u8; 64];
            keccak::shake256_2_into(&mut two, x, y);
            ok &= two.as_slice() == crypto_lib::hash::sha3::shake256(&joined, 64).as_slice();
            ok &= keccak::sha3_512_2(x, y).as_slice()
                == crypto_lib::hash::sha3::sha3_512(&joined).as_slice();
            println!("Keccak selftest: {}", if ok { "PASS" } else { "FAIL" });
            if !ok {
                std::process::exit(1);
            }
        }
        "bench" => {
            let block = vec![0xa5u8; 1024];
            let fast = timed(iters, || {
                std::hint::black_box(keccak::sha3_256(&block));
            });
            let reference = timed(iters, || {
                std::hint::black_box(crypto_lib::hash::sha3::sha3_256(&block));
            });
            println!("Keccak, SHA3-256 of 1024 bytes, {iters} iterations");
            println!("  pqc::fast:   {fast:>9.3} us");
            println!("  reference:   {reference:>9.3} us");
            println!("  speedup:     {:>9.2}x", reference / fast.max(1e-9));
        }
        other => {
            eprintln!("unknown op `{other}`; try demo, selftest or bench");
            std::process::exit(2);
        }
    }
}

/// A deterministic pseudo-random `F_p²` element, from a 64-bit seed.
fn fp2_from_seed(seed: u64) -> iso::Fp2 {
    let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    iso::Fp2::new(iso::Fp::from_u64(next()), iso::Fp::from_u64(next()))
}

/// A point of order exactly `2^e` on the supersingular curve `A = 0`, suitable
/// as the kernel generator of a `2^e`-isogeny.
///
/// `p + 1 = 3·2^324`, so multiplying any curve point by the cofactor 3 lands in
/// the 2-power torsion. `A = 0` (that is, `y² = x³ + x`) is supersingular
/// because `p ≡ 3 mod 4`.
///
/// Two conditions beyond "has 2-power order", both of which cost a retry when
/// they fail:
///
/// * **Full order.** `[2^{e-1}]K ≠ O` and `[2^e]K = O`.
/// * **Not above `(0,0)`.** `E[2^∞] ≅ Z/2^324 × Z/2^324` here, so a cyclic
///   subgroup of order `2^e` meets the 2-torsion in one of three subgroups, and
///   about a third of candidates sit above the 2-torsion point `(0,0)`. The
///   degree-4 isogeny formula is degenerate for those — which is exactly what
///   [`iso::isog4_kernel_ok`] tests — so they have to be discarded rather than
///   walked. Getting this wrong looks like "the chain returns `None`", which is
///   not obviously a *kernel choice* problem until you check.
fn two_power_point(curve: &iso::Curve, e: usize, seed: u64) -> Option<iso::PointX> {
    if !(2..=iso::TWO_TORSION_POWER).contains(&e) {
        return None;
    }
    for k in 0..128u64 {
        let x = fp2_from_seed(seed ^ (k << 32));
        let p = iso::PointX::from_affine(x);
        if !curve.is_on_curve(&p) {
            continue;
        }
        // Clear the cofactor 3, then drop to order exactly 2^e.
        let mut kp = iso::ladder_curve(&[3u64], 2, &p, curve);
        if iso::TWO_TORSION_POWER > e {
            kp = iso::xdbl_e_curve(&kp, iso::TWO_TORSION_POWER - e, curve);
        }
        if kp.is_infinity() {
            continue;
        }
        if iso::xdbl_e_curve(&kp, e - 1, curve).is_infinity() {
            continue; // order too small
        }
        if !iso::xdbl_e_curve(&kp, e, curve).is_infinity() {
            continue; // order too large, i.e. cofactor not cleared
        }
        // The first degree-4 step's kernel is this point doubled down to order
        // 4; it must not lie above (0,0).
        if e >= 2 {
            let order4 = iso::xdbl_e_curve(&kp, e - 2, curve);
            if !iso::isog4_kernel_ok(&order4, curve) {
                continue;
            }
        }
        return Some(kp);
    }
    None
}

fn isogeny(op: &str, e: usize, iters: u32) {
    let curve = iso::Curve::from_a(iso::Fp2::ZERO);
    match op {
        "info" => {
            println!("SQIsign-scale isogeny arithmetic (pqc::fast)");
            println!(
                "  p = 3·2^{} - 1, {} bits, {} limbs",
                iso::TWO_TORSION_POWER,
                iso::P_BITS,
                iso::NLIMBS
            );
            println!("  F_p² = F_p[i]/(i²+1); curves in x-only Montgomery (A24plus : C24) form");
            println!("  available 2-power torsion: 2^{}", iso::TWO_TORSION_POWER);
            println!(
                "  4-isogeny strategy costs: doubling {}, evaluation {}",
                iso::COST_DBL,
                iso::COST_EVAL
            );
            let a = curve.a_affine();
            println!(
                "  starting curve A = {:?} (supersingular: p ≡ 3 mod 4)",
                a.re.to_canonical()[0]
            );
        }
        "selftest" => {
            // Named stages, so a failure says which invariant broke rather
            // than only that something did.
            let mut stages: Vec<(&str, bool)> = Vec::new();

            let mut field = true;
            for s in 0..8u64 {
                let a = fp2_from_seed(s);
                let b = fp2_from_seed(s ^ 0xff);
                field &= a.mul(&b).sub(&b.mul(&a)).is_zero();
                field &= a.sqr().sub(&a.mul(&a)).is_zero();
                if !a.is_zero() {
                    field &= a.mul(&a.inv()).sub(&iso::Fp2::ONE).is_zero();
                }
                if let Some(r) = a.sqrt() {
                    field &= r.sqr().sub(&a).is_zero();
                }
            }
            stages.push(("F_p2 arithmetic identities", field));

            // xDBLADD must agree with xDBL and xADD done separately. A formula
            // identity, so it holds for any x-coordinates.
            let a24 = curve.normalised_a24();
            let mut ladder_ok = true;
            for s in 0..4u64 {
                let p = iso::PointX::from_affine(fp2_from_seed(s));
                let q = iso::PointX::from_affine(fp2_from_seed(s ^ 0xabc));
                let pmq = iso::PointX::from_affine(fp2_from_seed(s ^ 0xdef));
                let (d, a) = iso::xdblad(&p, &q, &pmq, &a24);
                ladder_ok &= d.same_x(&iso::xdbl_a24(&p, &a24));
                ladder_ok &= a.same_x(&iso::xadd(&p, &q, &pmq));
            }
            stages.push(("xDBLADD agrees with xDBL and xADD", ladder_ok));

            // A point of full 2-power order on the starting curve.
            let kernel = two_power_point(&curve, e.min(iso::TWO_TORSION_POWER), 1);
            stages.push(("point of full 2-power order found", kernel.is_some()));

            if let Some(k) = kernel {
                stages.push(("kernel generator lies on the curve", curve.is_on_curve(&k)));
                let e = e.min(iso::TWO_TORSION_POWER);
                let mut push = [iso::PointX::from_affine(fp2_from_seed(77))];
                let chain = iso::two_isogeny_chain(&curve, &k, e, &mut push);
                stages.push(("2^e isogeny chain completes", chain.is_some()));

                let n = e / 2;
                let strat = iso::optimal_strategy(n, iso::COST_DBL, iso::COST_EVAL);
                let mut p1 = [iso::PointX::from_affine(fp2_from_seed(77))];
                let mut p2 = p1;
                let c1 = iso::chain_4_strategy(&curve, &k, n, &mut p1, &strat);
                let c2 = iso::chain_4_naive(&curve, &k, n, &mut p2);
                match (c1, c2) {
                    (Some(a), Some(b)) => {
                        stages.push((
                            "strategy and naive walkers reach the same curve",
                            a.j_invariant().sub(&b.j_invariant()).is_zero(),
                        ));
                        stages.push(("both walkers push the same point", p1[0].same_x(&p2[0])));
                    }
                    _ => stages.push(("both 4-isogeny walkers complete", false)),
                }
            }

            let ok = stages.iter().all(|(_, v)| *v);
            for (name, v) in &stages {
                println!("  {:<48} {}", name, if *v { "ok" } else { "FAILED" });
            }
            println!("isogeny selftest: {}", if ok { "PASS" } else { "FAIL" });
            if !ok {
                std::process::exit(1);
            }
        }
        "strategy" => {
            let n = e.min(iso::TWO_TORSION_POWER) / 2;
            let strat = iso::optimal_strategy(n, iso::COST_DBL, iso::COST_EVAL);
            // The naive walker doubles down from the top each step: that is
            // n(n-1)/2 doublings plus n evaluations per pushed point.
            let naive = iso::COST_DBL * (n as u64 * (n as u64 - 1) / 2);
            let opt: u64 = strat.iter().map(|&s| iso::COST_DBL * s as u64).sum::<u64>()
                + iso::COST_EVAL * strat.len() as u64;
            println!(
                "optimal 4-isogeny strategy for {n} steps (2^{} chain)",
                n * 2
            );
            println!("  strategy length: {}", strat.len());
            println!("  naive walker cost:    {naive}");
            println!("  strategy cost:        {opt}");
            println!(
                "  saving:               {:.1}x",
                naive as f64 / opt.max(1) as f64
            );
        }
        "chain" => {
            let e = e.min(iso::TWO_TORSION_POWER);
            let Some(k) = two_power_point(&curve, e, 1) else {
                eprintln!("could not find a point of full 2-power order");
                std::process::exit(1);
            };
            let mut push = [iso::PointX::from_affine(fp2_from_seed(77))];
            let t = Instant::now();
            let codomain = iso::two_isogeny_chain(&curve, &k, e, &mut push);
            let dt = t.elapsed();
            match codomain {
                Some(c) => {
                    let j = c.j_invariant();
                    println!("2^{e}-isogeny walked in {:.3} ms", dt.as_secs_f64() * 1e3);
                    println!(
                        "  codomain j-invariant (re, low limb): {:?}",
                        j.re.to_canonical()[0]
                    );
                    println!("  pushed point still finite: {}", !push[0].is_infinity());
                }
                None => {
                    eprintln!("the chain hit a degenerate kernel");
                    std::process::exit(1);
                }
            }
        }
        "bench" => {
            let a = fp2_from_seed(1);
            let b = fp2_from_seed(2);
            let mul = timed(iters * 100, || {
                std::hint::black_box(a.mul(&b));
            });
            let sqr = timed(iters * 100, || {
                std::hint::black_box(a.sqr());
            });
            let inv = timed(iters, || {
                std::hint::black_box(a.inv());
            });
            let a24 = curve.normalised_a24();
            let p = iso::PointX::from_affine(fp2_from_seed(3));
            let dbl = timed(iters * 10, || {
                std::hint::black_box(iso::xdbl_a24(&p, &a24));
            });
            println!("isogeny arithmetic (pqc::fast), 326-bit p");
            println!("  Fp2 mul:   {mul:>9.4} us");
            println!("  Fp2 sqr:   {sqr:>9.4} us");
            println!(
                "  Fp2 inv:   {inv:>9.4} us  ({:.0}x a mul)",
                inv / mul.max(1e-12)
            );
            println!("  xDBL:      {dbl:>9.4} us");
            if let Some(k) = two_power_point(&curve, iso::TWO_TORSION_POWER, 1) {
                let n = iso::TWO_TORSION_POWER / 2;
                let strat = iso::optimal_strategy(n, iso::COST_DBL, iso::COST_EVAL);
                let mut push = [iso::PointX::from_affine(fp2_from_seed(77))];
                let t = Instant::now();
                let _ = iso::two_isogeny_chain_with_strategy(
                    &curve,
                    &k,
                    iso::TWO_TORSION_POWER,
                    &mut push,
                    &strat,
                );
                println!(
                    "  2^{} chain: {:>7.2} ms",
                    iso::TWO_TORSION_POWER,
                    t.elapsed().as_secs_f64() * 1e3
                );
            }
        }
        other => {
            eprintln!("unknown op `{other}`; try info, selftest, chain, strategy or bench");
            std::process::exit(2);
        }
    }
}
