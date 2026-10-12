//! # Root-finding decomposition on the NIST prime curves P-192 … P-521
//!
//! Runs the root-finding (direct-subtraction) point-decomposition oracle
//! on the five NIST prime curves — **P-192, P-224, P-256, P-384, P-521**,
//! with their real FIPS 186-4 parameters — and compares it head-to-head
//! against the Gröbner (`groebner_f4::buchberger`) decomposition on the
//! same ideal, over the same field, on the same targets.
//!
//! ## What is measured, per curve
//!
//! - **A. Head-to-head vs F4.** Fixed small factor base (default
//!   `|F| = 6`, the largest at which a Gröbner basis over a 192–521-bit
//!   field is still practical per target). For each target `R` the
//!   root-finder and `buchberger` on `{S₃(X, Y; x_R), ∏(X−x_i), ∏(Y−x_i)}`
//!   are both run and timed. Fail-closed cross-checks on every target:
//!   a found pair makes Semaev `S₃` vanish and makes *every* Gröbner
//!   basis element vanish; the two methods must agree on
//!   decomposable-vs-refuted.
//! - **B. Root-finder scaling.** `|F| = 64`, hundreds of random on-curve
//!   targets (the rejection path, which dominates in practice) plus
//!   positive controls `R = F_i ± F_j` re-verified in the group.
//! - **C. Per-curve specialisation (P-256).** The identical oracle on
//!   the crate's dedicated constant-time Montgomery field
//!   (`ecc::p256_field::P256FieldElement`) versus the generic
//!   arbitrary-precision field, to show what a per-curve field buys.
//!   P-192/224/384/521 run on the generic field here; their Solinas
//!   reductions are the same follow-on specialisation.
//!
//! Both root-finder paths use Montgomery's batched inversion (one field
//! inversion per target, amortised over the factor base).
//!
//! ```bash
//! cargo run --release --example semaev_rootfind_nist_bench
//! BENCH_OUT=research/semaev_rootfind_nist_20261012 \
//!   cargo run --release --example semaev_rootfind_nist_bench
//! # knobs: F4_FB=6 F4_CONSTRUCTED=3 F4_RANDOM=2 SCALE_FB=64 SCALE_RANDOM=200
//! ```
//!
//! ## This is not a security result
//!
//! Only the per-target *decomposition subroutine* is measured, on a
//! bounded synthetic factor base of small-`x` points. A fast decomposition
//! oracle says nothing about the ECDLP on these curves: the relation
//! collection and linear-algebra stages of index calculus remain
//! exponentially out of reach at these sizes, and the asymptotics are
//! unchanged. No key-recovery, sub-rho, or deployed-curve-security claim
//! is made or implied. Public parameters and synthetic known-answer
//! targets only; every decomposition is re-verified in the group.

use std::collections::HashMap;
use std::time::{Duration, Instant};

use num_bigint::{BigUint, RandBigInt};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crypto_lib::cryptanalysis::ec_index_calculus::{semaev_s3, sqrt_mod_p};
use crypto_lib::cryptanalysis::groebner_f4::{buchberger, reduce_basis, Ordering};
use crypto_lib::cryptanalysis::symmetrized_semaev::MPoly;
use crypto_lib::ecc::field::FieldElement;
use crypto_lib::ecc::p256_field::P256FieldElement;

// ── NIST P-curves (FIPS 186-4). All have a = −3. ───────────────────

struct Nist {
    name: &'static str,
    p_hex: &'static str,
    b_hex: &'static str,
}

const NIST: [Nist; 5] = [
    Nist {
        name: "P-192",
        p_hex: "fffffffffffffffffffffffffffffffeffffffffffffffff",
        b_hex: "64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1",
    },
    Nist {
        name: "P-224",
        p_hex: "ffffffffffffffffffffffffffffffff000000000000000000000001",
        b_hex: "b4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4",
    },
    Nist {
        name: "P-256",
        p_hex: "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
        b_hex: "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
    },
    Nist {
        name: "P-384",
        p_hex: "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff",
        b_hex: "b3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef",
    },
    Nist {
        name: "P-521",
        p_hex: "1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
        b_hex: "0051953eb9618e1c9a1f929a21a0b68540eea2da725b99b315f3b8b489918ef109e156193951ec7e937b1652c0bd3bb1bf073573df883d2c34f1ef451fd46b503f00",
    },
];

type Pt = (FieldElement, FieldElement);
type Rel = (usize, i8, usize, i8);

struct Curve {
    name: &'static str,
    p: BigUint,
    a: FieldElement,
    b: FieldElement,
}

impl Curve {
    fn new(n: &Nist) -> Self {
        let p = BigUint::parse_bytes(n.p_hex.as_bytes(), 16).unwrap();
        let a = fe(&(&p - BigUint::from(3u32)), &p);
        let b = fe(&BigUint::parse_bytes(n.b_hex.as_bytes(), 16).unwrap(), &p);
        Curve { name: n.name, p, a, b }
    }
    fn bits(&self) -> u64 {
        self.p.bits()
    }
    fn rhs(&self, x: &FieldElement) -> FieldElement {
        x.mul(x).mul(x).add(&self.a.mul(x)).add(&self.b)
    }
    /// Lift `x = 1, 2, …` to on-curve points until `size` are found.
    fn factor_base(&self, size: usize) -> (Vec<Pt>, HashMap<BigUint, usize>) {
        let mut fb = Vec::with_capacity(size);
        let mut xmap = HashMap::with_capacity(size);
        let mut x = 1u64;
        while fb.len() < size {
            let xfe = fe(&BigUint::from(x), &self.p);
            if let Some(y) = sqrt_mod_p(&self.rhs(&xfe).value, &self.p) {
                xmap.insert(xfe.value.clone(), fb.len());
                fb.push((xfe, fe(&y, &self.p)));
            }
            x += 1;
        }
        (fb, xmap)
    }
    fn random_point(&self, rng: &mut StdRng) -> Pt {
        loop {
            let xv = rng.gen_biguint_below(&self.p);
            if xv == BigUint::from(0u32) {
                continue;
            }
            let xfe = fe(&xv, &self.p);
            if let Some(y) = sqrt_mod_p(&self.rhs(&xfe).value, &self.p) {
                return (xfe, fe(&y, &self.p));
            }
        }
    }
}

fn fe(v: &BigUint, p: &BigUint) -> FieldElement {
    FieldElement::new(v.clone(), p.clone())
}

// ── Generic (arbitrary-precision) root-finding decomposition ────────

/// Montgomery's trick: invert every element of `ds` (all non-zero) with
/// a single field inversion and `3(k−1)` multiplications.
fn batch_inverse(ds: &[FieldElement]) -> Vec<FieldElement> {
    let k = ds.len();
    if k == 0 {
        return Vec::new();
    }
    let mut prefix = Vec::with_capacity(k);
    let mut acc = ds[0].clone();
    prefix.push(acc.clone());
    for d in &ds[1..] {
        acc = acc.mul(d);
        prefix.push(acc.clone());
    }
    let mut inv = acc.inv().expect("product of non-zero elements is invertible");
    let mut out = vec![FieldElement::zero(ds[0].modulus.clone()); k];
    for i in (1..k).rev() {
        out[i] = inv.mul(&prefix[i - 1]);
        inv = inv.mul(&ds[i]);
    }
    out[0] = inv;
    out
}

/// 2-decomposition `R = ε_i F_i + ε_j F_j` by direct subtraction: for
/// each `F_i`, compute `x(R ∓ F_i)` and look it up in the factor base.
/// This *is* the Semaev `S₃` root `x` of `S₃(x_R, x_i, x) = 0`, read
/// off from the group law instead of from a polynomial system.
fn big_decompose(fb: &[Pt], xmap: &HashMap<BigUint, usize>, xr: &FieldElement, yr: &FieldElement) -> Option<Rel> {
    let mut idx = Vec::with_capacity(fb.len());
    let mut ds = Vec::with_capacity(fb.len());
    for (i, (xi, _)) in fb.iter().enumerate() {
        let dx = xr.sub(xi);
        if !dx.is_zero() {
            idx.push(i);
            ds.push(dx);
        }
    }
    let invs = batch_inverse(&ds);
    for (k, &i) in idx.iter().enumerate() {
        let (xi, yi) = &fb[i];
        for (eps, num) in [(1i8, yr.add(yi)), (-1i8, yr.sub(yi))] {
            let lam = num.mul(&invs[k]);
            let x3 = lam.mul(&lam).sub(xr).sub(xi);
            if let Some(&j) = xmap.get(&x3.value) {
                if j >= i {
                    let y3 = lam.mul(&xr.sub(&x3)).sub(yr);
                    let tau = if y3 == fb[j].1 { 1i8 } else { -1 };
                    return Some((i, eps, j, tau));
                }
            }
        }
    }
    None
}

fn big_point_add(pt: &Pt, qt: &Pt, a: &FieldElement) -> Option<Pt> {
    let (x1, y1) = pt;
    let (x2, y2) = qt;
    if x1.value == x2.value {
        if y1.value == y2.value && !y1.is_zero() {
            let two = fe(&BigUint::from(2u32), &x1.modulus);
            let three = fe(&BigUint::from(3u32), &x1.modulus);
            let num = three.mul(&x1.mul(x1)).add(a);
            let den = two.mul(y1).inv()?;
            let lam = num.mul(&den);
            let x3 = lam.mul(&lam).sub(x1).sub(x2);
            let y3 = lam.mul(&x1.sub(&x3)).sub(y1);
            return Some((x3, y3));
        }
        return None;
    }
    let lam = y2.sub(y1).mul(&x2.sub(x1).inv()?);
    let x3 = lam.mul(&lam).sub(x1).sub(x2);
    let y3 = lam.mul(&x1.sub(&x3)).sub(y1);
    Some((x3, y3))
}

/// Fail-closed: re-add the claimed summands in the group and compare.
fn big_verify(fb: &[Pt], a: &FieldElement, r: &Pt, rel: Rel) -> bool {
    let neg = |pt: &Pt| (pt.0.clone(), pt.1.neg());
    let (i, ei, j, ej) = rel;
    let pi = if ei < 0 { neg(&fb[i]) } else { fb[i].clone() };
    let pj = if ej < 0 { neg(&fb[j]) } else { fb[j].clone() };
    match big_point_add(&pi, &pj, a) {
        Some(sum) => sum.0.value == r.0.value && sum.1.value == r.1.value,
        None => false,
    }
}

// ── P-256 specialised path (dedicated Montgomery field) ─────────────

type Pt256 = (P256FieldElement, P256FieldElement);

fn p256_batch_inverse(ds: &[P256FieldElement]) -> Vec<P256FieldElement> {
    let k = ds.len();
    if k == 0 {
        return Vec::new();
    }
    let mut prefix = Vec::with_capacity(k);
    let mut acc = ds[0];
    prefix.push(acc);
    for d in &ds[1..] {
        acc = acc.mul(d);
        prefix.push(acc);
    }
    let mut inv = acc.inv();
    let mut out = vec![P256FieldElement::ZERO; k];
    for i in (1..k).rev() {
        out[i] = inv.mul(&prefix[i - 1]);
        inv = inv.mul(&ds[i]);
    }
    out[0] = inv;
    out
}

fn p256_decompose(
    fb: &[Pt256],
    xmap: &HashMap<[u8; 32], usize>,
    xr: &P256FieldElement,
    yr: &P256FieldElement,
) -> Option<Rel> {
    let mut idx = Vec::with_capacity(fb.len());
    let mut ds = Vec::with_capacity(fb.len());
    for (i, (xi, _)) in fb.iter().enumerate() {
        let dx = xr.sub(xi);
        if !bool::from(dx.ct_is_zero()) {
            idx.push(i);
            ds.push(dx);
        }
    }
    let invs = p256_batch_inverse(&ds);
    for (k, &i) in idx.iter().enumerate() {
        let (xi, yi) = &fb[i];
        for (eps, num) in [(1i8, yr.add(yi)), (-1i8, yr.sub(yi))] {
            let lam = num.mul(&invs[k]);
            let x3 = lam.sqr().sub(xr).sub(xi);
            if let Some(&j) = xmap.get(&x3.to_bytes_be()) {
                if j >= i {
                    let y3 = lam.mul(&xr.sub(&x3)).sub(yr);
                    let tau = if bool::from(y3.ct_eq(&fb[j].1)) { 1i8 } else { -1 };
                    return Some((i, eps, j, tau));
                }
            }
        }
    }
    None
}

// ── Gröbner (F4) baseline on the same decomposition ideal ───────────

fn eval_mpoly(poly: &MPoly, point: &[FieldElement]) -> FieldElement {
    let mut acc = FieldElement::zero(poly.p.clone());
    for (exps, coef) in &poly.terms {
        let mut term = coef.clone();
        for (i, &e) in exps.iter().enumerate() {
            for _ in 0..e {
                term = term.mul(&point[i]);
            }
        }
        acc = acc.add(&term);
    }
    acc
}

/// Semaev `S₃(X, Y, x_R)` for `y² = x³ + a x + b` as a 2-variable `MPoly`:
/// `(X−Y)² x_R² − 2[(X+Y)(XY+a) + 2b] x_R + [(XY−a)² − 4b(X+Y)]`.
fn s3_in_xy(a: &FieldElement, b: &FieldElement, xr: &FieldElement) -> MPoly {
    let p = a.modulus.clone();
    let x = MPoly::variable(0, 2, p.clone());
    let y = MPoly::variable(1, 2, p.clone());
    let a_m = MPoly::constant(a.clone(), 2);
    let b_m = MPoly::constant(b.clone(), 2);
    let xr_m = MPoly::constant(xr.clone(), 2);
    let two = MPoly::constant(fe(&BigUint::from(2u32), &p), 2);
    let four = MPoly::constant(fe(&BigUint::from(4u32), &p), 2);
    let neg1 = FieldElement::one(p.clone()).neg();

    let diff = x.sub(&y);
    let term1 = diff.mul(&diff).mul(&xr_m).mul(&xr_m);
    let sum = x.add(&y);
    let prod = x.mul(&y);
    let inner = sum.mul(&prod.add(&a_m)).add(&two.mul(&b_m));
    let term2 = two.mul(&inner).scale(&neg1).mul(&xr_m);
    let part = prod.sub(&a_m);
    let term3 = part.mul(&part).sub(&four.mul(&b_m).mul(&sum));
    term1.add(&term2).add(&term3)
}

/// `∏_i (v_idx − x_i)`.
fn vanishing(fb: &[Pt], idx: usize, p: &BigUint) -> MPoly {
    let var = MPoly::variable(idx, 2, p.clone());
    let mut acc = MPoly::constant(FieldElement::one(p.clone()), 2);
    for (xi, _) in fb {
        acc = acc.mul(&var.sub(&MPoly::constant(xi.clone(), 2)));
    }
    acc
}

fn is_unit_ideal(reduced: &[MPoly]) -> bool {
    reduced.len() == 1
        && reduced[0].terms.len() == 1
        && reduced[0].terms.keys().next().map(|m| m.iter().all(|&e| e == 0)) == Some(true)
}

// ── helpers ─────────────────────────────────────────────────────────

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn fmt_dur(d: Duration) -> String {
    let us = d.as_secs_f64() * 1e6;
    if us < 1_000.0 {
        format!("{us:.1} µs")
    } else if us < 1_000_000.0 {
        format!("{:.2} ms", us / 1e3)
    } else {
        format!("{:.2} s", us / 1e6)
    }
}

fn fmt_speedup(slow: Duration, fast: Duration) -> String {
    if fast.as_secs_f64() <= 0.0 {
        return "—".into();
    }
    format!("{:.0}×", slow.as_secs_f64() / fast.as_secs_f64())
}

#[derive(Default)]
struct Json {
    rows: Vec<String>,
}

// ── main ────────────────────────────────────────────────────────────

fn main() {
    let f4_fb = env_usize("F4_FB", 6);
    let f4_constructed = env_usize("F4_CONSTRUCTED", 3);
    let f4_random = env_usize("F4_RANDOM", 2);
    let scale_fb = env_usize("SCALE_FB", 64);
    let scale_random = env_usize("SCALE_RANDOM", 200);
    let scale_controls = 100usize;

    println!("Root-finding decomposition on NIST P-192 / P-224 / P-256 / P-384 / P-521");
    println!("(decomposition subroutine only — NOT a security claim; see module docs)");
    println!(
        "knobs: F4_FB={f4_fb} F4_CONSTRUCTED={f4_constructed} F4_RANDOM={f4_random} SCALE_FB={scale_fb} SCALE_RANDOM={scale_random} controls={scale_controls}\n"
    );

    let mut json = Json::default();
    let mut all_beat_f4 = true;

    for n in NIST.iter() {
        let c = Curve::new(n);
        let mut rng = StdRng::seed_from_u64(0x5151_u64 ^ c.bits());
        println!("=== {} ({} -bit prime, a = −3, FIPS 186-4 b) ===", c.name, c.bits());

        // ── A. Head-to-head vs F4 ──────────────────────────────────
        let (fb, xmap) = c.factor_base(f4_fb);
        let vx = vanishing(&fb, 0, &c.p);
        let vy = vanishing(&fb, 1, &c.p);

        let mut targets: Vec<(Pt, bool)> = Vec::new();
        let mut tries = 0;
        while targets.iter().filter(|t| t.1).count() < f4_constructed && tries < 10_000 {
            tries += 1;
            let i = rng.gen_range(0..fb.len());
            let j = rng.gen_range(0..fb.len());
            if i == j {
                continue;
            }
            if let Some(r) = big_point_add(&fb[i], &fb[j], &c.a) {
                targets.push((r, true));
            }
        }
        for _ in 0..f4_random {
            targets.push((c.random_point(&mut rng), false));
        }

        let mut t_root = Duration::ZERO;
        let mut t_f4 = Duration::ZERO;
        let mut decomp = 0usize;
        let mut refuted = 0usize;
        for (r, constructed) in &targets {
            let t = Instant::now();
            let rel = big_decompose(&fb, &xmap, &r.0, &r.1);
            t_root += t.elapsed();

            let t = Instant::now();
            let s3 = s3_in_xy(&c.a, &c.b, &r.0);
            let basis = buchberger(&[s3, vx.clone(), vy.clone()], Ordering::Grevlex);
            let reduced = reduce_basis(&basis, Ordering::Grevlex);
            t_f4 += t.elapsed();
            let unit = is_unit_ideal(&reduced);

            match rel {
                Some(rel) => {
                    decomp += 1;
                    assert!(big_verify(&fb, &c.a, r, rel), "{}: relation fails in the group", c.name);
                    let (i, _, j, _) = rel;
                    let s = semaev_s3(&r.0, &fb[i].0, &fb[j].0, &c.a, &c.b);
                    assert!(s.is_zero(), "{}: S₃ must vanish on the found pair", c.name);
                    for poly in &reduced {
                        let v = eval_mpoly(poly, &[fb[i].0.clone(), fb[j].0.clone()]);
                        assert!(v.is_zero(), "{}: Gröbner element must vanish at the known root", c.name);
                    }
                    assert!(!unit, "{}: F4 said {{1}} but a decomposition exists", c.name);
                }
                None => {
                    assert!(!constructed, "{}: constructed sum must decompose", c.name);
                    assert!(unit, "{}: root-finder refuted but F4 did not return {{1}}", c.name);
                    refuted += 1;
                }
            }
        }
        let nt = targets.len() as u32;
        let per_root = t_root / nt;
        let per_f4 = t_f4 / nt;
        let beats = t_root < t_f4;
        all_beat_f4 &= beats;
        println!(
            "  A. vs F4   |F|={f4_fb:<3} targets={nt} (decomposed {decomp}, refuted {refuted}, all cross-checks pass)"
        );
        println!(
            "             root-find {:>10}/target   F4 {:>10}/target   speedup {}",
            fmt_dur(per_root),
            fmt_dur(per_f4),
            fmt_speedup(t_f4, t_root)
        );
        json.rows.push(format!(
            "{{\"curve\":\"{}\",\"bits\":{},\"section\":\"vs_f4\",\"fb\":{},\"targets\":{},\"decomposed\":{},\"refuted\":{},\"root_us_per_target\":{:.3},\"f4_us_per_target\":{:.3},\"speedup\":{:.1},\"all_checks_pass\":true}}",
            c.name, c.bits(), f4_fb, nt, decomp, refuted,
            per_root.as_secs_f64() * 1e6, per_f4.as_secs_f64() * 1e6,
            t_f4.as_secs_f64() / t_root.as_secs_f64().max(1e-12)
        ));

        // ── B. Root-finder scaling on a bigger factor base ────────
        let (fb, xmap) = c.factor_base(scale_fb);
        let randoms: Vec<Pt> = (0..scale_random).map(|_| c.random_point(&mut rng)).collect();
        let t = Instant::now();
        let mut hits = 0usize;
        for r in &randoms {
            if let Some(rel) = big_decompose(&fb, &xmap, &r.0, &r.1) {
                assert!(big_verify(&fb, &c.a, r, rel), "{}: random-target relation fails", c.name);
                hits += 1;
            }
        }
        let per_reject = t.elapsed() / scale_random as u32;

        let mut controls: Vec<Pt> = Vec::new();
        while controls.len() < scale_controls {
            let i = rng.gen_range(0..fb.len());
            let j = rng.gen_range(0..fb.len());
            if i == j {
                continue;
            }
            let neg_j = rng.gen_bool(0.5);
            let fj = if neg_j { (fb[j].0.clone(), fb[j].1.neg()) } else { fb[j].clone() };
            if let Some(r) = big_point_add(&fb[i], &fj, &c.a) {
                controls.push(r);
            }
        }
        // Time the decomposition only; verify (group re-add + S₃) afterwards.
        let t = Instant::now();
        let rels: Vec<Rel> = controls
            .iter()
            .map(|r| big_decompose(&fb, &xmap, &r.0, &r.1).expect("constructed sum must decompose"))
            .collect();
        let per_control = t.elapsed() / scale_controls as u32;
        let mut verified = 0usize;
        for (r, &rel) in controls.iter().zip(&rels) {
            assert!(big_verify(&fb, &c.a, r, rel), "{}: control relation fails", c.name);
            let (i, _, j, _) = rel;
            assert!(semaev_s3(&r.0, &fb[i].0, &fb[j].0, &c.a, &c.b).is_zero());
            verified += 1;
        }
        println!(
            "  B. scale   |F|={scale_fb:<3} random {scale_random} → {hits} hits, {:>10}/target (generic field); controls verified {verified}/{scale_controls}, {:>10}/target",
            fmt_dur(per_reject),
            fmt_dur(per_control)
        );
        json.rows.push(format!(
            "{{\"curve\":\"{}\",\"bits\":{},\"section\":\"scale_generic\",\"fb\":{},\"random_targets\":{},\"random_hits\":{},\"reject_us_per_target\":{:.3},\"controls_verified\":{},\"controls_total\":{},\"control_us_per_target\":{:.3}}}",
            c.name, c.bits(), scale_fb, scale_random, hits,
            per_reject.as_secs_f64() * 1e6, verified, scale_controls,
            per_control.as_secs_f64() * 1e6
        ));

        // ── C. Per-curve specialised field (P-256) ────────────────
        if c.name == "P-256" {
            let fb256: Vec<Pt256> = fb
                .iter()
                .map(|(x, y)| (P256FieldElement::from_biguint(&x.value), P256FieldElement::from_biguint(&y.value)))
                .collect();
            let xmap256: HashMap<[u8; 32], usize> =
                fb256.iter().enumerate().map(|(i, (x, _))| (x.to_bytes_be(), i)).collect();
            let to256 = |r: &Pt| (P256FieldElement::from_biguint(&r.0.value), P256FieldElement::from_biguint(&r.1.value));

            let randoms256: Vec<Pt256> = randoms.iter().map(to256).collect();
            let t = Instant::now();
            let mut hits256 = 0usize;
            for (k, r) in randoms256.iter().enumerate() {
                if let Some(rel) = p256_decompose(&fb256, &xmap256, &r.0, &r.1) {
                    assert!(big_verify(&fb, &c.a, &randoms[k], rel), "P-256 specialised relation fails");
                    hits256 += 1;
                }
            }
            let per_reject256 = t.elapsed() / scale_random as u32;
            assert_eq!(hits256, hits, "specialised and generic paths must agree on random targets");

            let controls256: Vec<Pt256> = controls.iter().map(to256).collect();
            let t = Instant::now();
            let rels256: Vec<Rel> = controls256
                .iter()
                .map(|r| p256_decompose(&fb256, &xmap256, &r.0, &r.1).expect("constructed sum must decompose"))
                .collect();
            let per_control256 = t.elapsed() / scale_controls as u32;
            let mut verified256 = 0usize;
            for (k, &rel) in rels256.iter().enumerate() {
                assert!(big_verify(&fb, &c.a, &controls[k], rel), "P-256 specialised control fails");
                verified256 += 1;
            }
            println!(
                "  C. special |F|={scale_fb:<3} P256FieldElement: random {:>10}/target ({} vs generic), controls {verified256}/{scale_controls} {:>10}/target ({} vs generic)",
                fmt_dur(per_reject256),
                fmt_speedup(per_reject, per_reject256),
                fmt_dur(per_control256),
                fmt_speedup(per_control, per_control256)
            );
            json.rows.push(format!(
                "{{\"curve\":\"P-256\",\"bits\":256,\"section\":\"scale_specialised\",\"field\":\"P256FieldElement\",\"fb\":{},\"random_targets\":{},\"random_hits\":{},\"reject_us_per_target\":{:.3},\"controls_verified\":{},\"controls_total\":{},\"control_us_per_target\":{:.3},\"speedup_vs_generic_reject\":{:.2}}}",
                scale_fb, scale_random, hits256, per_reject256.as_secs_f64() * 1e6,
                verified256, scale_controls, per_control256.as_secs_f64() * 1e6,
                per_reject.as_secs_f64() / per_reject256.as_secs_f64().max(1e-12)
            ));
        }
        println!();
    }

    println!(
        "Closing check — root-finding beats F4 on every NIST P-curve: {}",
        if all_beat_f4 { "YES" } else { "NO" }
    );
    println!("(per-target decomposition subroutine on a bounded factor base; not an ECDLP or security result)");
    assert!(all_beat_f4, "closing criterion not met");

    if let Ok(dir) = std::env::var("BENCH_OUT") {
        std::fs::create_dir_all(&dir).expect("create BENCH_OUT");
        let path = format!("{dir}/results.json");
        let body = format!("[\n  {}\n]\n", json.rows.join(",\n  "));
        std::fs::write(&path, body).expect("write results.json");
        println!("\nwrote {path}");
    }
}
