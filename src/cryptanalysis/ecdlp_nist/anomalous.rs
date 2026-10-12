//! Anomalous prime-field curves (`#E(F_p) = p`) and the
//! Smart–Semaev–Satoh–Araki attack that solves their ECDLP in polynomial
//! time.
//!
//! ## The attack
//!
//! Lift `E`, `P` and `Q = [k]P` to `Ẽ / Z_p` (any lift of the curve
//! equation that is not the canonical lift; `P̃`, `Q̃` by Hensel), and
//! compute `[p]P̃` and `[p]Q̃` in Jacobian coordinates modulo `p²` —
//! ring operations only, so precision is exact.  Because `#E(F_p) = p`,
//! both land in the formal group `Ê(pZ_p)`: their `Z` coordinate is
//! divisible by `p`.  The formal-group parameter `t = −x/y = −XZ/Y` is then
//! in `pZ_p`, the formal logarithm is `ψ(t) = t + O(t²) ≡ t (mod p²)`, and
//! `ψ` is a homomorphism `Ê(pZ_p) → pZ_p`, so
//!
//! ```text
//!     k ≡ ψ([p]Q̃) / ψ([p]P̃) ≡ (t_Q / p) · (t_P / p)⁻¹   (mod p).
//! ```
//!
//! With probability `1/p` the chosen lift of `E` is the canonical lift and
//! `t_P / p ≡ 0`; the attack then retries with another lift of `(a, b)`
//! (both coefficients moved by multiples of `p`: moving `b` alone leaves a
//! `j = 0` curve isomorphic to its canonical lift).  Every answer is verified by `[k]P = Q` on `E(F_p)` before it is
//! returned, so a wrong answer is impossible; a `None` is a genuine
//! failure (non-anomalous input, or `P` not of order `p`).
//!
//! ## Constructing anomalous curves
//!
//! `t = 1` means `4p = 1 + |D| v²` for the CM discriminant `D` of the
//! Frobenius.  For the class-number-one discriminants `D ∈ {−11, −19,
//! −43, −67, −163}` (`D ≡ 5 mod 8`, so that `(1 + |D|v²)/4` is odd; `D = −7`
//! gives an even value) the `j`-invariant is a known integer, so
//! [`generate`] picks `v` with `p = (1 + |D|v²)/4` prime at the requested
//! size, reduces the CM curve `y² = x³ + 3j(1728−j)x + 2j(1728−j)²`
//! modulo `p`, and takes whichever of it and its quadratic twist has
//! `[p]P = O` (the other has order `p + 2`).  This gives anomalous curves
//! of any size with no point counting, and is the fixture for the
//! "anomalous" branch of the solver: none of the NIST curves is anomalous
//! (their traces are far from `1`, see the audit), so the branch can only
//! be shown working on constructed instances.

use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::group::{EcdlpGroup, PrimeGroup};
use crate::cryptanalysis::curve_traits::arith::is_prime;
use crate::cryptanalysis::ec_index_calculus::sqrt_mod_p;
use crate::ecc::curve::CurveParams;
use crate::ecc::point::Point;

/// What one Smart-attack run did.
#[derive(Clone, Debug)]
pub struct SmartReport {
    /// The recovered `k`, verified by `[k]P = Q` on `E(F_p)`.
    pub scalar: Option<BigUint>,
    /// Lifts of `b` tried (`b + jp`, `j = 0, 1, …`) before one worked.
    pub lifts_tried: u32,
    pub elapsed_ms: u128,
    /// Why it failed, when it did.
    pub failure: Option<String>,
}

/// Jacobian point `(X : Y : Z)` over `Z / p²Z`.
#[derive(Clone, Debug)]
struct Jac {
    x: BigUint,
    y: BigUint,
    z: BigUint,
}

struct Zp2 {
    p: BigUint,
    p2: BigUint,
    a: BigUint,
}

impl Zp2 {
    #[inline]
    fn mul(&self, a: &BigUint, b: &BigUint) -> BigUint {
        a * b % &self.p2
    }
    #[inline]
    fn add(&self, a: &BigUint, b: &BigUint) -> BigUint {
        let s = a + b;
        if s >= self.p2 {
            s - &self.p2
        } else {
            s
        }
    }
    #[inline]
    fn sub(&self, a: &BigUint, b: &BigUint) -> BigUint {
        if a >= b {
            a - b
        } else {
            a + &self.p2 - b
        }
    }
    fn small(&self, k: u32, a: &BigUint) -> BigUint {
        a * BigUint::from(k) % &self.p2
    }

    /// `2P` (dbl-2007-bl with general `a`).
    fn double(&self, pt: &Jac) -> Jac {
        let xx = self.mul(&pt.x, &pt.x);
        let yy = self.mul(&pt.y, &pt.y);
        let yyyy = self.mul(&yy, &yy);
        let zz = self.mul(&pt.z, &pt.z);
        let xpyy = self.add(&pt.x, &yy);
        let s = self.small(2, &self.sub(&self.sub(&self.mul(&xpyy, &xpyy), &xx), &yyyy));
        let zz2 = self.mul(&zz, &zz);
        let m = self.add(&self.small(3, &xx), &self.mul(&self.a, &zz2));
        let t = self.sub(&self.mul(&m, &m), &self.small(2, &s));
        let y3 = self.sub(&self.mul(&m, &self.sub(&s, &t)), &self.small(8, &yyyy));
        let ypz = self.add(&pt.y, &pt.z);
        let z3 = self.sub(&self.sub(&self.mul(&ypz, &ypz), &yy), &zz);
        Jac { x: t, y: y3, z: z3 }
    }

    /// `P + Q` for `P ≠ ±Q` as `p`-adic points (add-2007-bl).
    fn add_pts(&self, p: &Jac, q: &Jac) -> Jac {
        let z1z1 = self.mul(&p.z, &p.z);
        let z2z2 = self.mul(&q.z, &q.z);
        let u1 = self.mul(&p.x, &z2z2);
        let u2 = self.mul(&q.x, &z1z1);
        let s1 = self.mul(&self.mul(&p.y, &q.z), &z2z2);
        let s2 = self.mul(&self.mul(&q.y, &p.z), &z1z1);
        let h = self.sub(&u2, &u1);
        let two_h = self.small(2, &h);
        let i = self.mul(&two_h, &two_h);
        let j = self.mul(&h, &i);
        let r = self.small(2, &self.sub(&s2, &s1));
        let v = self.mul(&u1, &i);
        let x3 = self.sub(&self.sub(&self.mul(&r, &r), &j), &self.small(2, &v));
        let y3 = self.sub(
            &self.mul(&r, &self.sub(&v, &x3)),
            &self.small(2, &self.mul(&s1, &j)),
        );
        let z1pz2 = self.add(&p.z, &q.z);
        let z3 = self.mul(
            &self.sub(&self.sub(&self.mul(&z1pz2, &z1pz2), &z1z1), &z2z2),
            &h,
        );
        Jac {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// `[k]P` by left-to-right double-and-add; the identity never occurs
    /// as an intermediate when `P` has order `p` on `E(F_p)` and `k ≤ p`.
    fn mul_pt(&self, pt: &Jac, k: &BigUint) -> Jac {
        let bits = k.bits();
        let mut acc = pt.clone();
        for i in (0..bits - 1).rev() {
            acc = self.double(&acc);
            if k.bit(i) {
                acc = self.add_pts(&acc, pt);
            }
        }
        acc
    }
}

fn inv_mod(a: &BigUint, m: &BigUint) -> Option<BigUint> {
    // Extended Euclid on BigUint via BigInt.
    use num_bigint::BigInt;
    use num_integer::Integer;
    let a = BigInt::from(a.clone());
    let m_i = BigInt::from(m.clone());
    let e = a.extended_gcd(&m_i);
    if !e.gcd.is_one() {
        return None;
    }
    let x = ((e.x % &m_i) + &m_i) % &m_i;
    x.to_biguint()
}

/// Hensel-lift `(x, y) ∈ E(F_p)` to `Ẽ : y² = x³ + ãx + b̃` modulo `p²`,
/// keeping `x̃ = x`.
fn hensel_lift(ctx: &Zp2, b: &BigUint, x: &BigUint, y: &BigUint) -> Option<(BigUint, BigUint)> {
    let p = &ctx.p;
    let p2 = &ctx.p2;
    let rhs = (x * x % p2 * x + &ctx.a * x + b) % p2;
    let y2 = y * y % p2;
    let f = (rhs + p2 - y2) % p2; // ≡ 0 mod p
    if (&f % p) != BigUint::zero() {
        return None;
    }
    let f_over_p = &f / p;
    let two_y_inv = inv_mod(&(BigUint::from(2u32) * y % p), p)?;
    let s = f_over_p * two_y_inv % p;
    let y_lift = (y + p * s) % p2;
    debug_assert_eq!(
        (&y_lift * &y_lift) % p2,
        (x * x % p2 * x + &ctx.a * x + b) % p2
    );
    Some((x.clone(), y_lift))
}

/// `t = −XZ/Y mod p²` of a Jacobian point in the formal group, returned as
/// `t / p mod p`; `None` if `Y` is not a unit or `Z` is not divisible by
/// `p` (the point did not land in the formal group: not anomalous).
fn formal_parameter_over_p(ctx: &Zp2, pt: &Jac) -> Result<BigUint, &'static str> {
    let p = &ctx.p;
    if (&pt.z % p) != BigUint::zero() {
        return Err(
            "[p]P̃ is not in the formal group: the curve is not anomalous or P has order ≠ p",
        );
    }
    let y_inv = inv_mod(&pt.y, &ctx.p2).ok_or("Y of [p]P̃ is not a unit mod p²")?;
    let xz = ctx.mul(&pt.x, &pt.z);
    let t = ctx.sub(&BigUint::zero(), &ctx.mul(&xz, &y_inv)); // −XZ/Y
    debug_assert_eq!(&t % p, BigUint::zero());
    Ok((t / p) % p)
}

/// Solve `Q = [k]P` on the anomalous curve `y² = x³ + ax + b` over `F_p`
/// (`#E(F_p) = p`).  Returns `k` verified on `E(F_p)`.
pub fn smart_attack(
    p: &BigUint,
    a: &BigUint,
    b: &BigUint,
    px: &BigUint,
    py: &BigUint,
    qx: &BigUint,
    qy: &BigUint,
) -> SmartReport {
    let t0 = std::time::Instant::now();
    let mut report = SmartReport {
        scalar: None,
        lifts_tried: 0,
        elapsed_ms: 0,
        failure: None,
    };
    let curve = CurveParams {
        name: "smart-target",
        p: p.clone(),
        a: a % p,
        b: b % p,
        gx: px % p,
        gy: py % p,
        n: p.clone(),
        h: 1,
    };
    let group = PrimeGroup::new(curve.clone());
    let p_pt = group.point(&curve.gx, &curve.gy);
    let q_pt = group.point(&(qx % p), &(qy % p));
    if !group.is_on_curve(&p_pt) || !group.is_on_curve(&q_pt) {
        report.failure = Some("P or Q is not on the curve".into());
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }
    if py.is_zero() || qy.is_zero() {
        report.failure = Some("P or Q has order 2 — not in a group of odd prime order".into());
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }
    let p2 = p * p;
    let verify = |k: &BigUint| group.mul(&p_pt, k) == q_pt;

    for j in 0..16u32 {
        report.lifts_tried = j + 1;
        // Lift j = 0 is (a, b) itself; later lifts move both coefficients
        // by pseudo-random multiples of p.  Moving b alone is not enough:
        // when a ≡ 0 (j-invariant 0) every (a, b + jp) is isomorphic to the
        // canonical lift, since the substitution u = 1 + mp sends b to
        // b(1 + 6mp) and fixes a.
        let (ra, rb) = if j == 0 {
            (BigUint::zero(), BigUint::zero())
        } else {
            let h = super::group::mix64(0x5A4D_u64.wrapping_mul(j as u64 + 1));
            (
                BigUint::from(h & 0xFFFF_FFFF) % p,
                BigUint::from(h >> 32) % p,
            )
        };
        let a_lift = (a % p) + p * ra;
        let b_lift = (b % p) + p * rb;
        let ctx = Zp2 {
            p: p.clone(),
            p2: p2.clone(),
            a: a_lift % &p2,
        };
        let Some((px_l, py_l)) = hensel_lift(&ctx, &b_lift, &curve.gx, &curve.gy) else {
            report.failure = Some("Hensel lift of P failed".into());
            break;
        };
        let Some((qx_l, qy_l)) = hensel_lift(&ctx, &b_lift, &(qx % p), &(qy % p)) else {
            report.failure = Some("Hensel lift of Q failed".into());
            break;
        };
        let one = BigUint::one();
        let pp = ctx.mul_pt(
            &Jac {
                x: px_l,
                y: py_l,
                z: one.clone(),
            },
            p,
        );
        let pq = ctx.mul_pt(
            &Jac {
                x: qx_l,
                y: qy_l,
                z: one,
            },
            p,
        );
        let u_p = match formal_parameter_over_p(&ctx, &pp) {
            Ok(u) => u,
            Err(e) => {
                report.failure = Some(e.into());
                break;
            }
        };
        let u_q = match formal_parameter_over_p(&ctx, &pq) {
            Ok(u) => u,
            Err(e) => {
                report.failure = Some(e.into());
                break;
            }
        };
        if u_p.is_zero() {
            // Canonical lift (probability 1/p): try another lift of b.
            continue;
        }
        let Some(inv) = inv_mod(&u_p, p) else {
            continue;
        };
        let k = u_q * inv % p;
        if verify(&k) {
            report.scalar = Some(k);
            break;
        }
        report.failure = Some("formal-log ratio failed verification".into());
        break;
    }
    if report.scalar.is_none() && report.failure.is_none() {
        report.failure = Some("every lift tried was canonical (astronomically unlikely)".into());
    }
    report.elapsed_ms = t0.elapsed().as_millis();
    report
}

/// Run the Smart attack on a [`PrimeGroup`] (its generator as `P`).
pub fn smart_attack_group(g: &PrimeGroup, target: &Point) -> SmartReport {
    let c = &g.curve;
    match target {
        Point::Infinity => SmartReport {
            scalar: Some(BigUint::zero()),
            lifts_tried: 0,
            elapsed_ms: 0,
            failure: None,
        },
        Point::Affine { x, y } => smart_attack(&c.p, &c.a, &c.b, &c.gx, &c.gy, &x.value, &y.value),
    }
}

/// Whether the prime-field curve has `#E(F_p) = p`, from its declared
/// `h · n`.
pub fn is_anomalous(curve: &CurveParams) -> bool {
    &curve.n * BigUint::from(curve.h) == curve.p
}

/// Class-number-one CM discriminants with `D ≡ 5 (mod 8)` (so that
/// `p = (1 + |D|v²)/4` is odd for odd `v`), and their integer
/// `j`-invariants (as `(|D|, |j|)`; every `j` here is negative).
const CM_TABLE: &[(u32, u128, &str)] = &[
    (11, 32768, "anomalous-cm-d11"),
    (19, 884_736, "anomalous-cm-d19"),
    (43, 884_736_000, "anomalous-cm-d43"),
    (67, 147_197_952_000, "anomalous-cm-d67"),
    (163, 262_537_412_640_768_000, "anomalous-cm-d163"),
];

/// An anomalous curve of about `bits` bits, deterministic in `seed`
/// (which also selects the CM discriminant).  `Err` only for `bits < 8`.
pub fn generate(bits: u32, seed: u64) -> Result<CurveParams, String> {
    if bits < 8 {
        return Err("anomalous curve: need bits ≥ 8".into());
    }
    let (d, j_abs, name) = CM_TABLE[(seed % CM_TABLE.len() as u64) as usize];
    let d_big = BigUint::from(d);
    // 4p ≈ 2^(bits+2) = 1 + D v²  ⇒  v ≈ √(2^(bits+2) / D).
    let target = (BigUint::one() << (bits + 2)) / &d_big;
    let mut v = target.sqrt();
    // Seed-dependent odd offset so different seeds give different primes.
    let span = (&v / 8u32).max(BigUint::one());
    v += BigUint::from(2 * (seed / CM_TABLE.len() as u64)) % (span * 2u32);
    if (&v % 2u32).is_zero() {
        v += BigUint::one();
    }
    let p = loop {
        let cand = (BigUint::one() + &d_big * &v * &v) / 4u32;
        if cand.bits() > bits as u64 + 2 {
            return Err("anomalous curve: no prime found near the requested size".into());
        }
        if is_prime(&cand) {
            break cand;
        }
        v += BigUint::from(2u32);
    };
    // j mod p for negative j.
    let j = (&p - BigUint::from(j_abs) % &p) % &p;
    let c = (BigUint::from(1728u32) + &p - &j) % &p; // 1728 − j
    let a = BigUint::from(3u32) * &j % &p * &c % &p;
    let b = BigUint::from(2u32) * &j % &p * &c % &p * &c % &p;

    let make = |a: &BigUint, b: &BigUint| -> Option<CurveParams> {
        let mut curve = CurveParams {
            name,
            p: p.clone(),
            a: a.clone(),
            b: b.clone(),
            gx: BigUint::zero(),
            gy: BigUint::zero(),
            n: p.clone(),
            h: 1,
        };
        let mut x = BigUint::from(seed % 7919 + 2);
        for _ in 0..4096 {
            let rhs = (&x * &x % &p * &x + a * &x + b) % &p;
            if let Some(y) = sqrt_mod_p(&rhs, &p) {
                if !y.is_zero() {
                    curve.gx = x.clone();
                    curve.gy = y;
                    let g = PrimeGroup::new(curve.clone());
                    let pt = g.generator();
                    if g.mul(&pt, &p) == Point::Infinity {
                        return Some(curve);
                    }
                    return None; // order p + 2: this is the wrong twist
                }
            }
            x += BigUint::one();
        }
        None
    };
    if let Some(c) = make(&a, &b) {
        return Ok(c);
    }
    // Quadratic twist by a non-residue: a' = a c², b' = b c³.
    let mut nr = BigUint::from(2u32);
    loop {
        if sqrt_mod_p(&nr, &p).is_none() {
            break;
        }
        nr += BigUint::one();
    }
    let a_t = &a * &nr % &p * &nr % &p;
    let b_t = &b * &nr % &p * &nr % &p * &nr % &p;
    make(&a_t, &b_t).ok_or_else(|| "anomalous curve: neither twist has order p (bug)".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::canonical_lift::find_anomalous_curve;
    use crate::cryptanalysis::ecdlp_nist::group::EcdlpGroup;

    fn plant_and_solve(curve: &CurveParams, k: &BigUint) -> SmartReport {
        let g = PrimeGroup::new(curve.clone());
        let q = g.mul(&g.generator(), k);
        smart_attack_group(&g, &q)
    }

    #[test]
    fn smart_attack_on_brute_forced_tiny_anomalous_curves() {
        let mut found = 0;
        for p in [
            11u64, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
        ] {
            let Some((a, b)) = find_anomalous_curve(p) else {
                continue;
            };
            found += 1;
            // Any affine point generates (prime order p).
            let pb = BigUint::from(p);
            let mut curve = CurveParams {
                name: "tiny-anomalous",
                p: pb.clone(),
                a: BigUint::from(a),
                b: BigUint::from(b),
                gx: BigUint::zero(),
                gy: BigUint::zero(),
                n: pb.clone(),
                h: 1,
            };
            let mut x = 0u64;
            loop {
                let rhs = BigUint::from((x * x * x + a * x + b) % p);
                if let Some(y) = sqrt_mod_p(&rhs, &pb) {
                    if !y.is_zero() {
                        curve.gx = BigUint::from(x);
                        curve.gy = y;
                        break;
                    }
                }
                x += 1;
            }
            assert!(is_anomalous(&curve));
            for k in 1..p {
                let rep = plant_and_solve(&curve, &BigUint::from(k));
                assert_eq!(rep.scalar, Some(BigUint::from(k)), "p={p} k={k}: {rep:?}");
            }
        }
        assert!(
            found >= 3,
            "expected several tiny anomalous curves, found {found}"
        );
    }

    #[test]
    fn cm_generated_curves_are_anomalous_and_break() {
        for (bits, seed) in [
            (40u32, 0u64),
            (64, 1),
            (96, 2),
            (128, 3),
            (160, 4),
            (256, 5),
        ] {
            let curve = generate(bits, seed).expect("construct");
            assert!(is_anomalous(&curve));
            assert!(curve.p.bits() >= bits as u64 - 1 && curve.p.bits() <= bits as u64 + 2);
            let g = PrimeGroup::new(curve.clone());
            assert!(g.is_on_curve(&g.generator()));
            let k = (&curve.p * 3u32 / 7u32) + BigUint::from(seed);
            let rep = plant_and_solve(&curve, &k);
            assert_eq!(rep.scalar, Some(k), "bits={bits} seed={seed}: {rep:?}");
        }
    }

    #[test]
    fn smart_attack_refuses_a_non_anomalous_curve() {
        let curve = CurveParams::secp256k1();
        let g = PrimeGroup::new(curve.clone());
        let q = g.mul(&g.generator(), &BigUint::from(5u32));
        let rep = smart_attack_group(&g, &q);
        assert!(rep.scalar.is_none());
        assert!(rep.failure.unwrap().contains("not anomalous"));
    }
}
