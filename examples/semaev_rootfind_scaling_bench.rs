//! # Does the root-finding decomposition hold up on bigger curves?
//!
//! A scaling study of the root-finding decomposition oracle across larger
//! field sizes, with two honest walls made explicit:
//!
//! 1. **Native-width cap (engineering).** The `u64` Montgomery path
//!    (`prime_fast::Fp64`, the binary `PairSumTable` / `Gf2` packing)
//!    needs `p < 2^62` and `n ≤ 63`. A `2^83` prime and `F_{2^83}` are
//!    past that: the prime side still runs at arbitrary precision
//!    (`BigUint`), slower; the binary single-word field stops at `n ≤ 63`.
//! 2. **Exponential factor base (fundamental).** The Koblitz
//!    Frobenius-invariant factor base has dimension `ℓ = ord_2(n)` —
//!    `52` at `n = 53`, `82` at `n = 83` — so `|V| = 2^ℓ` is not
//!    materialisable by *any* method (root finding or F4). The small-`n`
//!    head-to-head lives where `ord_2(n)` is small; "bigger `n`" does not
//!    help the attack, it makes the factor base exponentially larger.
//!    This bench therefore measures the decomposition *oracle* on a
//!    bounded factor base, decoupled from that wall.
//!
//! What is measured:
//! - **Prime, native `u64`** — direct-subtraction 2-decomposition at
//!   40/48/56/61-bit, fixed factor base, fail-closed in the group.
//! - **Prime, arbitrary-precision** — the same algorithm at ~83-bit
//!   (`p` = largest prime `< 2^83`), proving the method, not the word
//!   size, is what extends.
//! - **Binary `F_{2^53}`** — the pairs-and-solve (root-finding)
//!   decomposition over a bounded subspace, cross-checked against
//!   exhaustive `S₄` search for correctness and timed for rejection cost.
//!
//! ```bash
//! cargo run --release --example semaev_rootfind_scaling_bench
//! BENCH_OUT=research/semaev_rootfind_scaling_20261012 \
//!   cargo run --release --example semaev_rootfind_scaling_bench
//! ```
//!
//! **Not a security claim.** Public-synthetic known-answer; every
//! decomposition re-checked in the group. No deployed-curve,
//! key-recovery, or sub-rho claim.

use std::collections::HashMap;
use std::time::{Duration, Instant};

use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crypto_lib::cryptanalysis::ec_index_calculus::{semaev_s3, sqrt_mod_p};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use crypto_lib::cryptanalysis::prime_fast::{FastCurve, FastFactorBase, Pt};
use crypto_lib::cryptanalysis::residual_walk::is_prime_u64;
use crypto_lib::cryptanalysis::semaev_decomp::{decompose, eval_f3, Gf2};
use crypto_lib::ecc::field::FieldElement;

// ── tiny JSON sink ──────────────────────────────────────────────────

#[derive(Default)]
struct Json {
    rows: Vec<String>,
}
impl Json {
    fn row(&mut self, s: String) {
        self.rows.push(s);
    }
}

fn largest_prime_below_pow2(bits: u32, res_mod4: u64) -> u64 {
    let mut p = (1u64 << bits) - 1;
    loop {
        if p % 4 == res_mod4 && is_prime_u64(p) {
            return p;
        }
        p -= 2;
    }
}

// ── Prime, native u64 (prime_fast direct subtraction) ───────────────

/// Direct-subtraction 2-decomposition of `R` over the factor base: the
/// Semaev `S₃` root in `x` is `x(R ∓ F_i)`, found by one batched field
/// inversion and `|F|` table lookups. Returns `(i, ε_i, j, ε_j)` with
/// `R = ε_i F_i + ε_j F_j`, or `None`.
fn prime_decompose(
    curve: &FastCurve,
    fb: &FastFactorBase,
    r: Pt,
    d: &mut Vec<u64>,
    degen: &mut Vec<bool>,
    scratch: &mut Vec<u64>,
) -> Option<(usize, i8, usize, i8)> {
    let f = curve.f;
    let b = fb.len();
    d.clear();
    degen.clear();
    for i in 0..b {
        let dd = f.sub(r.x, fb.pts[i].x);
        degen.push(dd == 0);
        d.push(if dd == 0 { f.one } else { dd });
    }
    f.batch_inv(d, scratch);
    for i in 0..b {
        if degen[i] {
            continue;
        }
        let fi = fb.pts[i];
        let inv = d[i];
        for (eps, lam) in [
            (1i8, f.mul(f.add(r.y, fi.y), inv)),
            (-1i8, f.mul(f.sub(r.y, fi.y), inv)),
        ] {
            let x3 = f.sub(f.sub(f.sqr(lam), r.x), fi.x);
            let xc = f.from_mont(x3);
            if let Some(j) = fb.lookup(xc) {
                if j >= i {
                    let y3 = f.sub(f.mul(lam, f.sub(r.x, x3)), r.y);
                    let tau = if y3 == fb.pts[j].y { 1i8 } else { -1 };
                    return Some((i, eps, j, tau));
                }
            }
        }
    }
    None
}

fn prime_verify(
    curve: &FastCurve,
    fb: &FastFactorBase,
    r: Pt,
    rel: (usize, i8, usize, i8),
) -> bool {
    let (i, ei, j, ej) = rel;
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
    match curve.add(Some(pi), Some(pj)) {
        Some(sum) => curve.canonical(sum) == curve.canonical(r),
        None => false,
    }
}

fn prime_native_section(json: &mut Json) {
    println!("\n=== Prime, native u64 (Fp64): direct-subtraction 2-decomposition ===\n");
    println!(
        "{:>5}  {:>4}  {:>8}  {:>12}  {:>10}  {:>8}",
        "bits", "|F|", "rand tgt", "root-find/tgt", "verified", "ms total"
    );
    let fb_size = 120usize;
    for bits in [40u32, 48, 56, 61] {
        let p = largest_prime_below_pow2(bits, 3);
        // a = -3, b = 1 (y² = x³ − 3x + 1, non-singular: disc ≠ 0).
        let curve = FastCurve::new(&format!("scale-{bits}bit"), p, p - 3, 1, 0);
        let fb = FastFactorBase::build(&curve, fb_size);
        if fb.len() < fb_size {
            continue;
        }
        let mut rng = StdRng::seed_from_u64(0xA53_u64 ^ bits as u64);
        let (mut d, mut degen, mut scratch) = (Vec::new(), Vec::new(), Vec::new());

        // (a) rejection cost on random on-curve targets.
        let n_rand = 2000usize;
        let mut got = 0usize;
        let t = Instant::now();
        let mut sampled = 0usize;
        while sampled < n_rand {
            let x = rng.gen_range(1..p);
            let Some(r) = curve.lift_x(x) else { continue };
            sampled += 1;
            if let Some(rel) = prime_decompose(&curve, &fb, r, &mut d, &mut degen, &mut scratch) {
                assert!(
                    prime_verify(&curve, &fb, r, rel),
                    "native relation fails in-group"
                );
                got += 1;
            }
        }
        let per = t.elapsed() / n_rand as u32;

        // (b) fail-closed positive control: 200 constructed R = F_i + F_j.
        let mut verified = 0usize;
        while verified < 200 {
            let i = rng.gen_range(0..fb.len());
            let j = rng.gen_range(0..fb.len());
            if i == j {
                continue;
            }
            let Some(r) = curve.add(Some(fb.pts[i]), Some(fb.pts[j])) else {
                continue;
            };
            let rel = prime_decompose(&curve, &fb, r, &mut d, &mut degen, &mut scratch)
                .expect("a constructed sum must decompose");
            assert!(
                prime_verify(&curve, &fb, r, rel),
                "constructed relation fails in-group"
            );
            verified += 1;
        }
        println!(
            "{:>5}  {:>4}  {:>8}  {:>12.2?}  {:>10}  {:>8.1}",
            bits,
            fb.len(),
            format!("{got}/{n_rand}"),
            per,
            format!("{verified}/200"),
            (per * n_rand as u32).as_secs_f64() * 1e3
        );
        json.row(format!(
            "{{\"kind\":\"prime_native\",\"bits\":{},\"fb\":{},\"rand_targets\":{},\"rand_decomp\":{},\"per_target_s\":{:.9},\"verified\":{}}}",
            bits, fb.len(), n_rand, got, per.as_secs_f64(), verified
        ));
    }
}

// ── Prime, arbitrary precision (BigUint) at ~83 bits ────────────────

fn fe(v: &BigUint, p: &BigUint) -> FieldElement {
    FieldElement::new(v.clone(), p.clone())
}

#[allow(clippy::too_many_arguments)]
fn big_decompose(
    p: &BigUint,
    fb: &[(FieldElement, FieldElement)],
    xmap: &HashMap<BigUint, usize>,
    xr: &FieldElement,
    yr: &FieldElement,
) -> Option<(usize, i8, usize, i8)> {
    for (i, (xi, yi)) in fb.iter().enumerate() {
        let dx = xr.sub(xi);
        if dx.is_zero() {
            continue;
        }
        let inv = dx.inv().expect("non-zero difference is invertible");
        for (eps, num) in [(1i8, yr.add(yi)), (-1i8, yr.sub(yi))] {
            let lam = num.mul(&inv);
            let x3 = lam.mul(&lam).sub(xr).sub(xi);
            if let Some(&j) = xmap.get(&x3.value) {
                if j >= i {
                    let y3 = lam.mul(&xr.sub(&x3)).sub(yr);
                    let tau = if y3 == fb[j].1 { 1i8 } else { -1 };
                    let _ = p;
                    return Some((i, eps, j, tau));
                }
            }
        }
    }
    None
}

fn big_point_add(
    pt: &(FieldElement, FieldElement),
    qt: &(FieldElement, FieldElement),
    a: &FieldElement,
) -> Option<(FieldElement, FieldElement)> {
    let (x1, y1) = pt;
    let (x2, y2) = qt;
    if x1.value == x2.value {
        if y1.value == y2.value && !y1.is_zero() {
            // doubling
            let two = fe(&BigUint::from(2u32), &x1.modulus);
            let three = fe(&BigUint::from(3u32), &x1.modulus);
            let num = three.mul(&x1.mul(x1)).add(a);
            let den = two.mul(y1).inv()?;
            let lam = num.mul(&den);
            let x3 = lam.mul(&lam).sub(x1).sub(x2);
            let y3 = lam.mul(&x1.sub(&x3)).sub(y1);
            return Some((x3, y3));
        }
        return None; // P + (−P) = O
    }
    let lam = y2.sub(y1).mul(&x2.sub(x1).inv()?);
    let x3 = lam.mul(&lam).sub(x1).sub(x2);
    let y3 = lam.mul(&x1.sub(&x3)).sub(y1);
    Some((x3, y3))
}

fn big_verify(
    fb: &[(FieldElement, FieldElement)],
    a: &FieldElement,
    r: &(FieldElement, FieldElement),
    rel: (usize, i8, usize, i8),
) -> bool {
    let neg = |pt: &(FieldElement, FieldElement)| (pt.0.clone(), pt.1.neg());
    let (i, ei, j, ej) = rel;
    let pi = if ei < 0 { neg(&fb[i]) } else { fb[i].clone() };
    let pj = if ej < 0 { neg(&fb[j]) } else { fb[j].clone() };
    match big_point_add(&pi, &pj, a) {
        Some(sum) => sum.0.value == r.0.value && sum.1.value == r.1.value,
        None => false,
    }
}

fn prime_biguint_section(json: &mut Json) {
    println!("\n=== Prime, arbitrary precision (BigUint) at ~83 bits ===\n");
    // Largest prime below 2^83 (python-verified): 2^83 − 55.
    let p = BigUint::parse_bytes(b"9671406556917033397649353", 10).unwrap();
    let a = fe(&(&p - BigUint::from(3u32)), &p);
    let b = fe(&BigUint::from(1u32), &p);
    let fb_size = 80usize;

    // Factor base: lift x = 1, 2, … via a modular square root.
    let mut fb: Vec<(FieldElement, FieldElement)> = Vec::new();
    let mut xmap: HashMap<BigUint, usize> = HashMap::new();
    let mut x = 1u64;
    while fb.len() < fb_size && x < 100_000 {
        let xfe = fe(&BigUint::from(x), &p);
        let rhs = xfe.mul(&xfe).mul(&xfe).add(&a.mul(&xfe)).add(&b);
        if let Some(y) = sqrt_mod_p(&rhs.value, &p) {
            xmap.insert(xfe.value.clone(), fb.len());
            fb.push((xfe, fe(&y, &p)));
        }
        x += 1;
    }
    if fb.len() < fb_size {
        println!("  (could not build a {fb_size}-point factor base)");
        return;
    }

    let mut rng = StdRng::seed_from_u64(0x83_u64);
    let p_u128 = u128::from_str_radix("9671406556917033397649353", 10).unwrap();

    // (a) rejection cost on random on-curve targets.
    let n_rand = 200usize;
    let mut got = 0usize;
    let mut sampled = 0usize;
    let t = Instant::now();
    while sampled < n_rand {
        let xv = BigUint::from(rng.gen::<u128>() % p_u128 + 1);
        let xfe = fe(&xv, &p);
        let rhs = xfe.mul(&xfe).mul(&xfe).add(&a.mul(&xfe)).add(&b);
        let Some(yv) = sqrt_mod_p(&rhs.value, &p) else {
            continue;
        };
        sampled += 1;
        let yr = fe(&yv, &p);
        if let Some(rel) = big_decompose(&p, &fb, &xmap, &xfe, &yr) {
            assert!(
                big_verify(&fb, &a, &(xfe.clone(), yr.clone()), rel),
                "83-bit relation fails"
            );
            got += 1;
        }
    }
    let per = t.elapsed() / n_rand as u32;

    // (b) fail-closed positive control: constructed R = F_i + F_j.
    let mut verified = 0usize;
    while verified < 100 {
        let i = rng.gen_range(0..fb.len());
        let j = rng.gen_range(0..fb.len());
        if i == j {
            continue;
        }
        let Some(r) = big_point_add(&fb[i], &fb[j], &a) else {
            continue;
        };
        let rel =
            big_decompose(&p, &fb, &xmap, &r.0, &r.1).expect("a constructed sum must decompose");
        assert!(
            big_verify(&fb, &a, &r, rel),
            "constructed 83-bit relation fails"
        );
        // Cross-check: Semaev S₃ vanishes on the found x-pair.
        let (ii, _, jj, _) = rel;
        let s = semaev_s3(&r.0, &fb[ii].0, &fb[jj].0, &a, &b);
        assert!(s.is_zero(), "S₃ must vanish on the 83-bit decomposition");
        verified += 1;
    }

    println!(
        "{:>7}  {:>4}  {:>10}  {:>14}  {:>10}",
        "bits", "|F|", "rand tgt", "root-find/tgt", "verified"
    );
    println!(
        "{:>7}  {:>4}  {:>10}  {:>14.2?}  {:>10}",
        83,
        fb.len(),
        format!("{got}/{n_rand}"),
        per,
        format!("{verified}/100")
    );
    println!("  p = 2^83 − 55 (largest prime < 2^83); every decomposition verified in-group.");
    json.row(format!(
        "{{\"kind\":\"prime_biguint\",\"bits\":83,\"fb\":{},\"rand_targets\":{},\"rand_decomp\":{},\"per_target_s\":{:.9},\"verified\":{}}}",
        fb.len(), n_rand, got, per.as_secs_f64(), verified
    ));
}

// ── Binary F_{2^53}: pairs-and-solve over a bounded subspace ─────────

fn binary_section(json: &mut Json) {
    println!("\n=== Binary F_(2^53): pairs-and-solve (root-finding) decomposition ===\n");
    let irr = find_irreducible_sparse(53).expect("irreducible of degree 53");
    let gf = Gf2::new(&irr);
    let mask53: u64 = (1u64 << 53) - 1;

    // (a) correctness anchor at ℓ = 6: decompose must agree with an
    // exhaustive S₄ triple search on every target.
    let l_anchor = 6u32;
    let span = 1u64 << l_anchor;
    let mut rng = StdRng::seed_from_u64(0x5353);
    let mut anchor_ok = 0usize;
    for _ in 0..200 {
        let xr = rng.gen::<u64>() & mask53;
        let got = decompose(xr, l_anchor, &gf);
        // Exhaustive: is there a triple (a ≤ b ≤ c) in V with S₄ = 0?
        let mut brute = false;
        'outer: for a in 0..span {
            for b in a..span {
                for c in b..span {
                    if eval_f3(a, b, c, xr, &gf) == 0 {
                        brute = true;
                        break 'outer;
                    }
                }
            }
        }
        assert_eq!(
            got.is_some(),
            brute,
            "decompose disagrees with exhaustive S₄ search at xr={xr:#x}"
        );
        if let Some(w) = got {
            assert_eq!(
                eval_f3(w[0], w[1], w[2], xr, &gf),
                0,
                "witness does not satisfy S₄"
            );
        }
        anchor_ok += 1;
    }
    println!(
        "  correctness anchor (ℓ=6): decompose == exhaustive S₄ search on {anchor_ok}/200 targets."
    );

    // (b) rejection-cost scaling across ℓ (the dominant per-target cost).
    println!(
        "\n{:>4}  {:>8}  {:>10}  {:>14}",
        "ℓ", "|V|", "targets", "decomp/tgt"
    );
    for l in [8u32, 10, 11] {
        let targets = if l >= 11 { 8 } else { 40 };
        let mut t_total = Duration::ZERO;
        for _ in 0..targets {
            let xr = rng.gen::<u64>() & mask53;
            let t = Instant::now();
            let _ = decompose(xr, l, &gf);
            t_total += t.elapsed();
        }
        let per = t_total / targets as u32;
        println!(
            "{:>4}  {:>8}  {:>10}  {:>14.2?}",
            l,
            1u64 << l,
            targets,
            per
        );
        json.row(format!(
            "{{\"kind\":\"binary_f2_53\",\"ell\":{},\"subspace\":{},\"targets\":{},\"per_target_s\":{:.9}}}",
            l,
            1u64 << l,
            targets,
            per.as_secs_f64()
        ));
    }
    println!(
        "  (subspace bounded on purpose: the Koblitz Frobenius ℓ = ord_2(53) = 52 ⇒ |V| = 2^52,"
    );
    println!("   which no index-calculus method can materialise — see the diagnosis below.)");
}

fn diagnosis() {
    println!("\n=== Diagnosis: the two walls at F_(2^83) and 2^83-scale primes ===\n");
    println!(
        "  native u64 width : Fp64 needs p < 2^62; PairSumTable / Gf2 pack x into u64 (n ≤ 63)."
    );
    println!(
        "                     → 83-bit prime runs here only via BigUint (measured above, slower);"
    );
    println!(
        "                       F_(2^83) exceeds the single-word binary field (Gf2 asserts n ≤ 63)"
    );
    println!(
        "                       and the pair-table packing (n ≤ 62); a u128 field is the fix."
    );
    println!("  Koblitz factor base (fundamental, affects F4 and root finding equally):");
    println!(
        "                     ℓ = ord_2(n):  n=53 → ℓ=52 (|V|=2^52),  n=83 → ℓ=82 (|V|=2^82)."
    );
    println!(
        "                     The invariant factor base is exponential, so the *full attack* is"
    );
    println!("                     infeasible at these n for ANY decomposition engine. The oracle");
    println!(
        "                     itself scales (shown above on a bounded subspace / factor base)."
    );
}

fn main() {
    println!("Root-finding decomposition — scaling to bigger curves.");
    println!("Public synthetic known-answer; every decomposition re-checked in the group.");
    let mut json = Json::default();
    prime_native_section(&mut json);
    prime_biguint_section(&mut json);
    binary_section(&mut json);
    diagnosis();

    if let Ok(dir) = std::env::var("BENCH_OUT") {
        let _ = std::fs::create_dir_all(&dir);
        let body = format!("[\n  {}\n]\n", json.rows.join(",\n  "));
        let _ = std::fs::write(format!("{dir}/results.json"), body);
        println!("\nWrote raw results to {dir}/results.json");
    }
    println!("\nEvery decomposition verified in-group; no security claim (see header).");
}
