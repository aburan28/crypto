//! Cross-checks for the isogeny evaluation module.
//!
//! - Field, polynomial and ladder arithmetic against naive references.
//! - √élu against Vélu (codomain and images) on CM instances with a fully
//!   rational `ℓ`-torsion (`j = −3375`, `p = (t² + 7v²)/4`, `ℓ | v`,
//!   `t ≡ 2 mod ℓ`, as in round-4 lead 16.3), for `ℓ ∈ {3, 5, 7, 11, 101,
//!   1009}`; codomain `j` against the library's independent
//!   `isogeny_walk::kernel::velu_codomain` (and `verify_kernel`, and
//!   `Φ_ℓ(j, j') = 0` for small `ℓ`); images have the right order on the
//!   codomain and the x-only map respects differential addition.
//! - `ℓ ∈ {2, 3}`: images on the codomain, homomorphism, certified dual with
//!   `φ̂∘φ = [ℓ]` and `φ∘φ̂ = [ℓ]`, the library codomain for `ℓ = 3`,
//!   `Φ_ℓ(j, j') = 0`, and `j ∈ {0, 1728}` domains.
//! - Counters: positive, deterministic, exact per-image Vélu cost, and the
//!   uncounted field gives the same results.

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use super::counted_fp::{is_prime_u64, NoCount, OpCounter, OpTally, PrimeField};
use super::montgomery::{batch_to_affine, is_kernel_generator, ladder, xadd, MontgomeryCurve, XZ};
use super::poly::{self, MonicDivisor, ProductTree};
use super::small_degree::{isomorphisms, kernel_x_coordinates, ShortWeierstrass, SmallIsogeny};
use super::sqrt_velu::sqrt_velu;
use super::velu::{kernel_polynomial, velu};
use super::IsogenyEvalError;
use crate::cryptanalysis::isogeny_walk::curve::Model;
use crate::cryptanalysis::isogeny_walk::field::Field as LibField;
use crate::cryptanalysis::isogeny_walk::{kernel as lib_kernel, modpoly, poly as lib_poly};

fn random_prime(rng: &mut StdRng, bits: u32) -> u64 {
    loop {
        let p = rng.gen_range(1u64 << (bits - 1)..1u64 << bits) | 1;
        if is_prime_u64(p) {
            return p;
        }
    }
}

fn lib_field(p: u64) -> LibField {
    LibField::new(&BigUint::from(p)).expect("prime")
}

fn lib_u64(f: &LibField, v: &crate::cryptanalysis::isogeny_walk::field::Fe) -> u64 {
    f.to_big(v).to_u64().expect("fits")
}

/// `Φ_ℓ(j1, j2) mod p` from the library's modular polynomial.
fn phi_at(p: u64, ell: u64, j1: u64, j2: u64) -> u64 {
    let lf = lib_field(p);
    let phi = modpoly::modular_polynomial(&lf, ell).expect("Φ_ℓ");
    let py = phi.at_x(&lf, &lf.from_u64(j1));
    lib_u64(&lf, &lib_poly::eval(&lf, &py, &lf.from_u64(j2)))
}

// --------------------------------------------------------------- field ----

#[test]
fn field_rejects_bad_moduli_and_inverts() {
    assert!(PrimeField::new(3).is_none());
    assert!(PrimeField::new(91).is_none());
    assert!(PrimeField::new((1u64 << 63) + 29).is_none());
    assert!(PrimeField::new(2_305_843_009_213_693_951).is_some()); // 2^61 − 1
    assert_eq!(std::mem::size_of::<NoCount>(), 0);
    let mut rng = StdRng::seed_from_u64(1);
    for bits in [10, 40, 61, 63] {
        let f = PrimeField::new(random_prime(&mut rng, bits)).unwrap();
        assert_eq!(f.inv(0), None);
        for _ in 0..50 {
            let a = rng.gen_range(1..f.p());
            assert_eq!(f.mul(a, f.inv(a).unwrap()), 1);
            let s = f.sqr(a);
            let r = f.sqrt(s).unwrap();
            assert_eq!(f.sqr(r), s);
            assert_eq!(f.pow(a, (f.p() - 1) as u128), 1);
        }
    }
}

// ---------------------------------------------------------- polynomials ----

fn naive_mul<C: OpCounter>(f: &PrimeField<C>, a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            out[i + j] = f.add(out[i + j], f.mul(x, y));
        }
    }
    out
}

#[test]
fn karatsuba_matches_schoolbook() {
    let mut rng = StdRng::seed_from_u64(2);
    for bits in [20, 62, 63] {
        let f = PrimeField::new(random_prime(&mut rng, bits)).unwrap();
        for (la, lb) in [
            (1, 1),
            (5, 30),
            (25, 25),
            (60, 61),
            (300, 40),
            (97, 200),
            (128, 128),
        ] {
            let a: Vec<u64> = (0..la).map(|_| rng.gen_range(0..f.p())).collect();
            let b: Vec<u64> = (0..lb).map(|_| rng.gen_range(0..f.p())).collect();
            assert_eq!(poly::mul(&f, &a, &b), naive_mul(&f, &a, &b), "{la}x{lb}");
        }
    }
}

#[test]
fn remainder_and_multipoint_evaluation() {
    let mut rng = StdRng::seed_from_u64(3);
    let f = PrimeField::new(random_prime(&mut rng, 62)).unwrap();
    for (dg, db) in [(10, 3), (100, 40), (300, 64), (80, 79)] {
        let mut b: Vec<u64> = (0..db).map(|_| rng.gen_range(0..f.p())).collect();
        b.push(1);
        let q: Vec<u64> = (0..=dg - db).map(|_| rng.gen_range(0..f.p())).collect();
        let r: Vec<u64> = (0..db).map(|_| rng.gen_range(0..f.p())).collect();
        let a = poly::add(&f, &poly::mul(&f, &q, &b), &r);
        let mut div = MonicDivisor::new(b).unwrap();
        let mut want = r.clone();
        poly::trim(&mut want);
        assert_eq!(div.rem(&f, &a), want);
        assert_eq!(div.rem(&f, &a), want, "cached reciprocal");
    }
    for n in [1usize, 2, 7, 37, 100] {
        let xs: Vec<u64> = (0..n).map(|_| rng.gen_range(0..f.p())).collect();
        let g: Vec<u64> = (0..150).map(|_| rng.gen_range(0..f.p())).collect();
        let mut tree = ProductTree::new(&f, &xs);
        let got = tree.evaluate(&f, &g);
        let want: Vec<u64> = xs.iter().map(|&x| poly::eval(&f, &g, x)).collect();
        assert_eq!(got, want, "n = {n}");
        let root = tree.root();
        assert!(xs.iter().all(|&x| poly::eval(&f, &root, x) == 0));
    }
    assert!(ProductTree::new(&f, &[]).evaluate(&f, &[1, 2]).is_empty());
}

// ------------------------------------------------- CM instances (λ = 1) ----

/// A Montgomery curve with `E[ℓ] ⊂ E(F_p)`, its order `n`, its twist's
/// order, and the x-coordinate of a point of order `ℓ`.
struct Inst {
    p: u64,
    curve: MontgomeryCurve,
    n: u128,
    n_twist: u128,
    kernel_x: u64,
}

fn valuation(mut n: u128, l: u128) -> (u32, u128) {
    let mut e = 0;
    while n.is_multiple_of(l) {
        n /= l;
        e += 1;
    }
    (e, n)
}

/// `j = −3375` (CM by `Z[(1+√−7)/2]`) over `p = (t² + 7v²)/4` with `ℓ | v`
/// and `t ≡ 2 (mod ℓ)`, `p ∈ [2^(bits−1), 2^bits)`.
fn cm_instance(ell: u64, bits: u32, rng: &mut StdRng) -> Inst {
    let lo: u128 = 1 << (bits + 1);
    let hi: u128 = 1 << (bits + 2);
    let l = ell as i128;
    'outer: loop {
        let wmax = (((hi - 1) / 7) as f64).sqrt() as u128 / ell as u128;
        let w = rng.gen_range(1..=wmax.max(1)) as u64;
        let v = ell * w;
        let sv = 7u128 * (v as u128) * (v as u128);
        if sv >= hi {
            continue;
        }
        let tmax = ((hi - 1 - sv) as f64).sqrt() as i128;
        let tmin = if sv >= lo {
            0
        } else {
            ((lo - sv) as f64).sqrt().ceil() as i128
        };
        let r = if (2 % l) % 2 == (v % 2) as i128 {
            2 % l
        } else {
            2 % l + l
        };
        let span = (tmax - tmin) / (2 * l);
        if span < 1 {
            continue;
        }
        let k = rng.gen_range(0..=span);
        let mut t = tmin - tmin.rem_euclid(2 * l) + r + 2 * l * k;
        if t < tmin {
            t += 2 * l;
        }
        let four_p = (t * t) as u128 + sv;
        if t > tmax || !four_p.is_multiple_of(4) || four_p < lo || four_p >= hi {
            continue;
        }
        let p = (four_p / 4) as u64;
        if !is_prime_u64(p) || p < 1 << (bits - 1) {
            continue;
        }
        let f = PrimeField::new(p).unwrap();
        let n = (p as i128 + 1 - t) as u128;
        let n_twist = (p as i128 + 1 + t) as u128;
        // Roots of x³ − 35x − 98 = (x − 7)(x² + 7x + 14).
        let mut roots = vec![7u64];
        if let Some(sq) = f.sqrt(f.from_i64(-7)) {
            let inv2 = f.inv(2).unwrap();
            roots.push(f.mul(f.sub(sq, 7), inv2));
            roots.push(f.mul(f.sub(f.neg(sq), 7), inv2));
        }
        for alpha in roots {
            let w = f.add(f.mul_small(3, f.sqr(alpha)), f.from_i64(-35));
            let Some(s) = f.sqrt(w).filter(|&s| s != 0) else {
                continue;
            };
            let a = f.div(f.mul_small(3, alpha), s).unwrap();
            let Ok(curve) = MontgomeryCurve::new(&f, a) else {
                continue;
            };
            assert_eq!(curve.j_invariant(&f), f.from_i64(-3375));
            let (e, m) = valuation(n, ell as u128);
            assert!(e >= 2, "ℓ² | N for λ = 1");
            for _ in 0..64 {
                let x = rng.gen_range(2..p);
                let mut q = ladder(&f, &XZ::affine(x), m, &curve.a24);
                if q.is_infinity() {
                    continue;
                }
                for _ in 0..e {
                    let r = ladder(&f, &q, ell as u128, &curve.a24);
                    if r.is_infinity() {
                        let kernel_x = q.to_affine(&f).unwrap();
                        return Inst {
                            p,
                            curve,
                            n,
                            n_twist,
                            kernel_x,
                        };
                    }
                    q = r;
                }
            }
            continue 'outer;
        }
    }
}

/// The order (`n` or the twist's) that kills `x` on `curve`.
fn side(f: &PrimeField, curve: &MontgomeryCurve, x: u64, inst: &Inst) -> u128 {
    if ladder(f, &XZ::affine(x), inst.n, &curve.a24).is_infinity() {
        inst.n
    } else {
        assert!(ladder(f, &XZ::affine(x), inst.n_twist, &curve.a24).is_infinity());
        inst.n_twist
    }
}

/// The library's codomain `j` for the kernel `⟨P⟩` (via the `B = 1` short
/// Weierstrass model), and the library kernel certificate for small `ℓ`.
fn library_codomain_j(f: &PrimeField, inst: &Inst, ell: u64) -> u64 {
    let (aw, bw, shift) = inst.curve.short_weierstrass(f);
    let d = ((ell - 1) / 2) as usize;
    let mults: Vec<XZ<u64>> = (1..=d as u128)
        .map(|i| ladder(f, &XZ::affine(inst.kernel_x), i, &inst.curve.a24))
        .collect();
    let xs = batch_to_affine(f, &mults).unwrap();
    let lf = lib_field(inst.p);
    let mut h = vec![lf.one()];
    for x in xs {
        let root = lf.from_u64(f.add(x, shift));
        h = lib_poly::mul(&lf, &h, &vec![lf.neg(&root), lf.one()]);
    }
    let e = Model {
        a: lf.from_u64(aw),
        b: lf.from_u64(bw),
    };
    let cod = lib_kernel::velu_codomain(&lf, &e, &h);
    if ell <= 11 {
        let certified = lib_kernel::verify_kernel(&lf, &e, ell, &h).expect("library certifies");
        assert_eq!(certified, cod);
    }
    lib_u64(&lf, &cod.j(&lf).expect("nonsingular"))
}

#[test]
fn sqrt_velu_agrees_with_velu_and_library() {
    let mut rng = StdRng::seed_from_u64(163);
    for (k, &ell) in [3u64, 5, 7, 11, 101, 1009].iter().enumerate() {
        for bits in [50u32, if k % 2 == 0 { 61 } else { 63 }] {
            let inst = cm_instance(ell, bits, &mut rng);
            let f = PrimeField::new(inst.p).unwrap();
            assert!(is_kernel_generator(&f, &inst.curve, ell, inst.kernel_x));
            // Random points on either side, a kernel point, and x = 0.
            let mut xs: Vec<u64> = (0..4).map(|_| rng.gen_range(2..inst.p)).collect();
            xs.push(inst.kernel_x);
            xs.push(0);
            let slow = velu(&f, &inst.curve, ell, inst.kernel_x, &xs).unwrap();
            let (fast, shape) = sqrt_velu(&f, &inst.curve, ell, inst.kernel_x, &xs).unwrap();
            assert_eq!(slow, fast, "ℓ = {ell}, p = {}", inst.p);
            if ell >= 5 {
                assert_eq!(
                    shape.i_len * shape.j_len * 2 + shape.k_len,
                    ((ell - 1) / 2) as usize
                );
            }
            assert_eq!(slow.images[4], None, "kernel point maps to O");
            assert_eq!(slow.images[5], Some(0), "(0, 0) maps to x = 0");
            let cod = slow.codomain;
            let j_new = cod.j_invariant(&f);
            assert_eq!(j_new, library_codomain_j(&f, &inst, ell), "ℓ = {ell}");
            if ell <= 7 {
                assert_eq!(phi_at(inst.p, ell, inst.curve.j_invariant(&f), j_new), 0);
            }
            // Images lie on the codomain's Kummer line with the same order.
            for (&x, img) in xs[..4].iter().zip(&slow.images) {
                let n = side(&f, &inst.curve, x, &inst);
                let xi = img.expect("non-kernel point");
                assert!(ladder(&f, &XZ::affine(xi), n, &cod.a24).is_infinity());
            }
        }
    }
}

#[test]
fn x_only_images_respect_differential_addition() {
    let mut rng = StdRng::seed_from_u64(5);
    for ell in [5u64, 11, 101] {
        let inst = cm_instance(ell, 55, &mut rng);
        let f = PrimeField::new(inst.p).unwrap();
        let c = &inst.curve;
        let g = XZ::affine(rng.gen_range(2..inst.p));
        let (a, b) = (rng.gen_range(1000..1u128 << 40), rng.gen_range(1..1000u128));
        let pts = [a, b, a + b, a - b].map(|k| ladder(&f, &g, k, &c.a24).to_affine(&f).unwrap());
        for use_sqrt in [false, true] {
            let out = if use_sqrt {
                sqrt_velu(&f, c, ell, inst.kernel_x, &pts).unwrap().0
            } else {
                velu(&f, c, ell, inst.kernel_x, &pts).unwrap()
            };
            let im: Vec<XZ<u64>> = out.images.iter().map(|x| XZ::affine(x.unwrap())).collect();
            let sum = xadd(&f, &im[0], &im[1], &im[3]);
            assert!(sum.same_as(&f, &im[2]), "ℓ = {ell}");
        }
    }
}

#[test]
fn invalid_inputs_are_errors_not_panics() {
    let mut rng = StdRng::seed_from_u64(6);
    let inst = cm_instance(7, 40, &mut rng);
    let f = PrimeField::new(inst.p).unwrap();
    assert_eq!(
        velu(&f, &inst.curve, 9, inst.kernel_x, &[]),
        Err(IsogenyEvalError::InvalidDegree(9))
    );
    assert_eq!(
        sqrt_velu(&f, &inst.curve, 2, inst.kernel_x, &[]).map(|r| r.1),
        Err(IsogenyEvalError::InvalidDegree(2))
    );
    assert!(MontgomeryCurve::new(&f, 2).is_err());
    assert!(MontgomeryCurve::new(&f, f.p() - 2).is_err());
    // A 7-torsion point is not a 101-kernel: no panic either way.
    assert!(!is_kernel_generator(&f, &inst.curve, 101, inst.kernel_x));
    let _ = sqrt_velu(&f, &inst.curve, 101, inst.kernel_x, &[5]);
    let _ = velu(&f, &inst.curve, 101, inst.kernel_x, &[5]);
    let h = kernel_polynomial(&f, &inst.curve, 7, inst.kernel_x).unwrap();
    assert_eq!(poly::degree(&h), Some(3));
    assert_eq!(poly::eval(&f, &h, inst.kernel_x), 0);
}

// ------------------------------------------------------- counters ----

#[test]
fn counters_are_positive_deterministic_and_free_when_off() {
    let mut rng = StdRng::seed_from_u64(7);
    let inst = cm_instance(101, 61, &mut rng);
    let xs = [12345u64, 678910];
    let plain = PrimeField::new(inst.p).unwrap();
    let counted = PrimeField::counted(inst.p).unwrap();
    let want = sqrt_velu(&plain, &inst.curve, 101, inst.kernel_x, &xs).unwrap();
    let mut tallies = Vec::new();
    for _ in 0..2 {
        counted.reset_counts();
        let got = sqrt_velu(&counted, &inst.curve, 101, inst.kernel_x, &xs).unwrap();
        assert_eq!(got, want);
        tallies.push(counted.tally());
    }
    assert_eq!(tallies[0], tallies[1]);
    assert!(tallies[0].mul > 0 && tallies[0].inv > 0);
    // Vélu: each extra image costs exactly 4d + 4 M and one I.
    let d = 50u64;
    let start = counted.tally();
    velu(&counted, &inst.curve, 101, inst.kernel_x, &xs[..1]).unwrap();
    let one = counted.tally().since(start);
    let start = counted.tally();
    velu(&counted, &inst.curve, 101, inst.kernel_x, &xs).unwrap();
    let two = counted.tally().since(start);
    assert_eq!(
        two.since(one),
        OpTally {
            mul: 4 * d + 4,
            inv: 1
        }
    );
    assert_eq!(two.mul_equivalents(91), two.mul + 91 * two.inv);
}

// ------------------------------------------------- ℓ ∈ {2, 3} ----

fn random_point<C: OpCounter>(
    f: &PrimeField<C>,
    e: &ShortWeierstrass,
    rng: &mut StdRng,
) -> (u64, u64) {
    e.point_from_x(f, rng.gen_range(0..f.p()), 1 << 20).unwrap()
}

fn check_small(f: &PrimeField, e: &ShortWeierstrass, ell: u64, x0: u64, rng: &mut StdRng) {
    let phi = SmallIsogeny::new(f, e, ell, x0).unwrap();
    let cod = *phi.codomain();
    let dual = phi.dual(f).unwrap();
    assert_eq!(dual.domain(), &cod);
    assert_eq!(dual.codomain(), e, "dual lands on the same model");
    assert_eq!(dual.degree(), ell);
    assert_eq!(phi_at(f.p(), ell, e.j_invariant(f), cod.j_invariant(f)), 0);
    if ell == 3 {
        let lf = lib_field(f.p());
        let lib = lib_kernel::velu_codomain(
            &lf,
            &Model {
                a: lf.from_u64(e.a),
                b: lf.from_u64(e.b),
            },
            &vec![lf.from_u64(f.neg(x0)), lf.one()],
        );
        assert_eq!((lib_u64(&lf, &lib.a), lib_u64(&lf, &lib.b)), (cod.a, cod.b));
    }
    for _ in 0..6 {
        let p = Some(random_point(f, e, rng));
        let r = Some(random_point(f, e, rng));
        let (ip, ir) = (phi.eval(f, p), phi.eval(f, r));
        assert!(cod.contains(f, ip) && cod.contains(f, ir));
        assert_eq!(
            phi.eval(f, e.add(f, p, r)),
            cod.add(f, ip, ir),
            "homomorphism"
        );
        assert_eq!(dual.eval(f, ip), e.mul(f, p, ell as u128), "dual∘φ = [ℓ]");
        let s = Some(random_point(f, &cod, rng));
        assert_eq!(
            phi.eval(f, dual.eval(f, s)),
            cod.mul(f, s, ell as u128),
            "φ∘dual = [ℓ]"
        );
    }
    // Kernel points map to O.
    if let Some(y) = f.sqrt(e.rhs(f, x0)) {
        assert_eq!(phi.eval(f, Some((x0, y))), None);
    }
}

#[test]
fn small_degree_isogenies_and_certified_duals() {
    let mut rng = StdRng::seed_from_u64(165);
    for bits in [30u32, 62] {
        let p = random_prime(&mut rng, bits);
        let f = PrimeField::new(p).unwrap();
        for ell in [2u64, 3] {
            let mut done = 0;
            while done < 3 {
                let (a, b) = (rng.gen_range(1..p), rng.gen_range(1..p));
                let Ok(e) = ShortWeierstrass::new(&f, a, b) else {
                    continue;
                };
                for x0 in kernel_x_coordinates(&f, &e, ell).unwrap() {
                    check_small(&f, &e, ell, x0, &mut rng);
                    done += 1;
                }
            }
        }
    }
}

#[test]
fn small_degree_special_j_and_errors() {
    let mut rng = StdRng::seed_from_u64(9);
    // p ≡ 1 (mod 12): j = 0 and j = 1728 have extra automorphisms.
    let p = loop {
        let p = random_prime(&mut rng, 40);
        if p % 12 == 1 {
            break p;
        }
    };
    let f = PrimeField::new(p).unwrap();
    // j = 1728: y² = x³ + ax, kernel (0, 0).
    let e1728 = ShortWeierstrass::new(&f, 5, 0).unwrap();
    assert_eq!(e1728.j_invariant(&f), 1728);
    check_small(&f, &e1728, 2, 0, &mut rng);
    // j = 0: y² = x³ + b with b = −c³ has the rational 2-torsion x = c.
    let c = 11u64;
    let e0 = ShortWeierstrass::new(&f, 0, f.neg(f.mul(f.sqr(c), c))).unwrap();
    assert_eq!(e0.j_invariant(&f), 0);
    for x0 in kernel_x_coordinates(&f, &e0, 2).unwrap() {
        check_small(&f, &e0, 2, x0, &mut rng);
    }
    for x0 in kernel_x_coordinates(&f, &e0, 3).unwrap() {
        check_small(&f, &e0, 3, x0, &mut rng);
    }
    assert!(isomorphisms(&f, &e0, &e0).len() >= 2);
    assert_eq!(isomorphisms(&f, &e0, &e1728), Vec::<u64>::new());
    // Errors.
    assert_eq!(
        SmallIsogeny::new(&f, &e1728, 5, 0),
        Err(IsogenyEvalError::InvalidDegree(5))
    );
    assert_eq!(
        SmallIsogeny::new(&f, &e1728, 2, 1),
        Err(IsogenyEvalError::NotAKernelPoint)
    );
    assert!(ShortWeierstrass::new(&f, 0, 0).is_err());
}

#[test]
fn small_degree_counts_are_deterministic() {
    let p = 1_000_000_007u64;
    let f = PrimeField::counted(p).unwrap();
    let e = ShortWeierstrass::new(&f, 0, f.neg(8)).unwrap(); // x = 2 is 2-torsion
    let phi = SmallIsogeny::new(&f, &e, 2, 2).unwrap();
    let pt = e.point_from_x(&f, 3, 1000);
    let start = f.tally();
    let img = phi.eval(&f, pt);
    assert_eq!(f.tally().since(start), OpTally { mul: 7, inv: 1 });
    assert!(phi.codomain().contains(&f, img));
}
