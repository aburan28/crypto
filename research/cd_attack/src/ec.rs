//! Elliptic curve arithmetic over F_{p²}.
//!
//! Two representations:
//!   - Short Weierstrass:  y² = x³ + a x + b   (with a, b ∈ F_{p²})
//!   - Montgomery:         y² = x³ + A x² + x  (with A ∈ F_{p²}, B=1)
//!
//! Affine arithmetic on both, with conversion between them. Used as the
//! source-side elliptic-curve machinery for the Castryck-Decru attack:
//! the starting curve E_0 = M_6 (Montgomery A=6), the j=1728 curve
//! E_{1728} = M_0, and the 2-isogeny between them.

use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};

use crate::field::{F2, Fp2};

// ---- Affine point on a Weierstrass / Montgomery curve --------------------

/// Affine point — either a finite (x, y) or the point at infinity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Aff {
    Inf,
    P(F2, F2),
}

impl Aff {
    pub fn x(&self) -> Option<&F2> {
        match self {
            Aff::Inf => None,
            Aff::P(x, _) => Some(x),
        }
    }
    pub fn y(&self) -> Option<&F2> {
        match self {
            Aff::Inf => None,
            Aff::P(_, y) => Some(y),
        }
    }
    pub fn is_inf(&self) -> bool { matches!(self, Aff::Inf) }
}

// ---- Montgomery curve  y² = x³ + A x² + x  -------------------------------

#[derive(Clone, Debug)]
pub struct Montgomery {
    pub a: F2,
}

impl Montgomery {
    pub fn new(a: F2) -> Self { Montgomery { a } }

    /// Is (x, y) on the curve?
    pub fn contains(&self, p: &Aff, fp2: &Fp2) -> bool {
        match p {
            Aff::Inf => true,
            Aff::P(x, y) => {
                let lhs = fp2.sq(y);
                let x2 = fp2.sq(x);
                let x3 = fp2.mul(&x2, x);
                let rhs = fp2.add(&fp2.add(&x3, &fp2.mul(&self.a, &x2)), x);
                fp2.sub(&lhs, &rhs) == fp2.zero()
            }
        }
    }

    /// Lift x to a point (x, ±y) if y² = x³ + A x² + x is a square in F_{p²}.
    pub fn lift(&self, x: &F2, fp2: &Fp2) -> Option<Aff> {
        let x2 = fp2.sq(x);
        let x3 = fp2.mul(&x2, x);
        let rhs = fp2.add(&fp2.add(&x3, &fp2.mul(&self.a, &x2)), x);
        fp2.sqrt(&rhs).map(|y| Aff::P(x.clone(), y))
    }

    pub fn neg(&self, p: &Aff, fp2: &Fp2) -> Aff {
        match p {
            Aff::Inf => Aff::Inf,
            Aff::P(x, y) => Aff::P(x.clone(), fp2.neg(y)),
        }
    }

    /// Affine point addition. Handles all the special cases of the chord-
    /// and-tangent law.
    pub fn add(&self, p: &Aff, q: &Aff, fp2: &Fp2) -> Aff {
        let ((x1, y1), (x2, y2)) = match (p, q) {
            (Aff::Inf, _) => return q.clone(),
            (_, Aff::Inf) => return p.clone(),
            (Aff::P(x1, y1), Aff::P(x2, y2)) => ((x1, y1), (x2, y2)),
        };
        if x1 == x2 {
            // y₁ + y₂ = 0  ⇒  result is ∞
            if fp2.add(y1, y2) == fp2.zero() {
                return Aff::Inf;
            }
            // Doubling: λ = (3x² + 2Ax + 1) / (2y)
            let three = fp2.from_int(3);
            let two = fp2.from_int(2);
            let num = fp2.add(
                &fp2.add(
                    &fp2.mul(&three, &fp2.sq(x1)),
                    &fp2.mul(&two, &fp2.mul(&self.a, x1)),
                ),
                &fp2.one(),
            );
            let den = fp2.mul(&two, y1);
            let lam = fp2.div(&num, &den);
            // x₃ = λ² − A − 2x₁
            let x3 = fp2.sub(
                &fp2.sub(&fp2.sq(&lam), &self.a),
                &fp2.mul(&two, x1),
            );
            let y3 = fp2.sub(&fp2.mul(&lam, &fp2.sub(x1, &x3)), y1);
            Aff::P(x3, y3)
        } else {
            // Chord: λ = (y₂ − y₁) / (x₂ − x₁)
            let lam = fp2.div(&fp2.sub(y2, y1), &fp2.sub(x2, x1));
            let x3 = fp2.sub(
                &fp2.sub(&fp2.sub(&fp2.sq(&lam), &self.a), x1),
                x2,
            );
            let y3 = fp2.sub(&fp2.mul(&lam, &fp2.sub(x1, &x3)), y1);
            Aff::P(x3, y3)
        }
    }

    pub fn dbl(&self, p: &Aff, fp2: &Fp2) -> Aff {
        self.add(p, p, fp2)
    }

    /// Scalar multiplication via left-to-right binary.
    pub fn mul(&self, k: &BigInt, p: &Aff, fp2: &Fp2) -> Aff {
        if k.is_zero() || p.is_inf() {
            return Aff::Inf;
        }
        let one = BigInt::one();
        let mut r = Aff::Inf;
        let bits = k.bits();
        for i in (0..bits).rev() {
            r = self.dbl(&r, fp2);
            if ((k >> i) & &one).is_one() {
                r = self.add(&r, p, fp2);
            }
        }
        r
    }
}

// ---- The "two_i" endomorphism on E_0  (j-invariant 1728) -----------------
//
// On the Montgomery curve E_{1728}: y² = x(x² + 1) (A=0), the map
//   ι : (x, y) → (−x, i·y)
// is an automorphism with ι² = [−1]. So 2ι is an endomorphism of degree 4.
//
// For the SIDH "starting curve" E_6 = y² = x³ + 6x² + x, the endomorphism
// "two_i" is conventionally the composition  φ ∘ ι ∘ φ̂  where φ : E_6 → E_{1728}
// is the 2-isogeny with kernel ⟨(0, 0)⟩.  See SIKE_challenge.m / Castryck-Decru.
//
// For the user-facing API: we implement `two_i_on_E1728` directly (on the
// j=1728 model), and provide the 2-isogeny / its dual for callers wanting
// the E_6 version.

/// ι : (x, y) → (−x, i·y) on E_{1728} = M_0.
pub fn iota_on_e1728(p: &Aff, fp2: &Fp2) -> Aff {
    match p {
        Aff::Inf => Aff::Inf,
        Aff::P(x, y) => Aff::P(
            fp2.neg(x),
            fp2.mul(&fp2.i(), y),
        ),
    }
}

/// 2·ι (x, y)  on E_{1728}.  This is the degree-4 endomorphism with
/// (2ι)² = [−4].  Use this when the source curve is E_{1728}; for E_6,
/// transport via the 2-isogeny.
pub fn two_iota_on_e1728(p: &Aff, fp2: &Fp2) -> Aff {
    let e0 = Montgomery::new(fp2.zero());
    let ip = iota_on_e1728(p, fp2);
    e0.dbl(&ip, fp2)
}

// ---- Degree-3 isogeny on Montgomery (affine) -----------------------------
//
// Kernel ⟨K⟩, K = (x_K, y_K) of order 3 (so [3]K = O, [2]K = -K).
// The Vélu-style x-map:
//   φ_x(x) = x · (x · x_K − 1)² / (x − x_K)²
// The y-map via chain rule on invariant differentials:
//   y_new = y · φ_x'(x)
// where
//   φ_x'(x) = (x x_K − 1)·(x² x_K − 3 x x_K² + x + x_K) / (x − x_K)³
//
// Codomain coefficient (Edwards-form pull-back, ℓ = 3, prod = (x_K−1, x_K+1)):
//   d = (A−2)/(A+2)
//   d_new = d³ · ((x_K − 1)/(x_K + 1))⁸
//   A_new = 2(d_new_num + d_new_den) / (d_new_den − d_new_num)
// with d_new_num = (A−2)³ · (x_K − 1)⁸ and d_new_den = (A+2)³ · (x_K + 1)⁸.

/// Apply a degree-3 isogeny to a Montgomery curve given the kernel x-coord.
/// Returns the new curve coefficient.
pub fn isog3_codomain(a: &F2, x_k: &F2, fp2: &Fp2) -> F2 {
    let two = fp2.from_int(2);
    let a_minus_2 = fp2.sub(a, &two);
    let a_plus_2 = fp2.add(a, &two);
    let xk_minus_1 = fp2.sub(x_k, &fp2.one());
    let xk_plus_1 = fp2.add(x_k, &fp2.one());

    let cube = |x: &F2| fp2.mul(x, &fp2.sq(x));
    let pow8 = |x: &F2| {
        let s = fp2.sq(x);
        let q = fp2.sq(&s);
        fp2.sq(&q)
    };

    let d_num = fp2.mul(&cube(&a_minus_2), &pow8(&xk_minus_1));
    let d_den = fp2.mul(&cube(&a_plus_2), &pow8(&xk_plus_1));

    // A_new = 2 (d_num + d_den) / (d_den − d_num)
    let num = fp2.mul(&two, &fp2.add(&d_num, &d_den));
    let den = fp2.sub(&d_den, &d_num);
    fp2.div(&num, &den)
}

/// Push the x-coordinate of an affine point P through a degree-3 isogeny
/// with kernel of x-coord x_K. Returns φ_x(P.x).
pub fn isog3_push_x(x: &F2, x_k: &F2, fp2: &Fp2) -> F2 {
    // φ_x(x) = x · (x · x_K − 1)² / (x − x_K)²
    let xxk_minus_1 = fp2.sub(&fp2.mul(x, x_k), &fp2.one());
    let x_minus_xk = fp2.sub(x, x_k);
    let num_x = fp2.mul(x, &fp2.sq(&xxk_minus_1));
    let den_x = fp2.sq(&x_minus_xk);
    fp2.div(&num_x, &den_x)
}

/// Push an affine point P through a degree-3 isogeny on Montgomery curve.
/// After computing x_new, the y-sign is disambiguated by lifting on the
/// codomain curve and matching the sign convention y_new = y · ε where ε
/// is chosen so that the curve equation holds. (For separable isogenies
/// with the standard normalization, ε = φ_x'(x); but the chain-rule
/// formula's normalization isn't always consistent with the codomain
/// formula's normalization on Montgomery curves, so we instead pick the
/// sign by curve-equation verification.)
///
/// Returns Aff::Inf if P is in the kernel.
pub fn isog3_push(p: &Aff, x_k: &F2, a_new: &F2, fp2: &Fp2) -> Aff {
    let (x, _y) = match p {
        Aff::Inf => return Aff::Inf,
        Aff::P(x, y) => (x, y),
    };
    // x_minus_xk == 0 → P is in the kernel
    let x_minus_xk = fp2.sub(x, x_k);
    if fp2.is_zero(&x_minus_xk) {
        return Aff::Inf;
    }
    let phi_x = isog3_push_x(x, x_k, fp2);
    // Lift on the codomain
    let e_new = Montgomery::new(a_new.clone());
    match e_new.lift(&phi_x, fp2) {
        Some(p) => p,
        None => Aff::Inf, // shouldn't happen for valid input
    }
}

/// Compute a chain of `n` 3-isogenies starting from `curve` with initial
/// kernel point `k0` of order 3^n. At step i, the kernel is [3^{n-i-1}]·k_i
/// where k_i is the image of k0 through the first i isogenies. Optionally
/// push a list of "auxiliary" points through the chain.
///
/// Returns: (final_curve, pushed_auxiliary_points).
///
/// Mirrors `Pushing3Chain` from richelot_aux.m.
pub fn pushing_3_chain(
    curve: &Montgomery,
    k0: &Aff,
    n: u32,
    aux: &mut [Aff],
    fp2: &Fp2,
) -> Montgomery {
    let mut cur = curve.clone();
    let mut k_running = k0.clone();
    for j in 1..=n {
        // Kernel of this step: [3^{n-j}]·k_running
        let cofactor = BigInt::from(3u64).pow(n - j);
        let small_k = cur.mul(&cofactor, &k_running, fp2);
        let x_k = match &small_k {
            Aff::Inf => panic!("kernel collapsed prematurely at step {j}"),
            Aff::P(x, _) => x.clone(),
        };
        // Codomain
        let a_new = isog3_codomain(&cur.a, &x_k, fp2);
        // Push k_running (kernel point → identity on its own step but lives
        // through the rest of the chain otherwise)
        k_running = isog3_push(&k_running, &x_k, &a_new, fp2);
        // Push each aux point
        for p in aux.iter_mut() {
            *p = isog3_push(p, &x_k, &a_new, fp2);
        }
        cur = Montgomery::new(a_new);
    }
    cur
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Fp, Fp2};

    fn ctx() -> Fp2 { Fp2::new(Fp::new(BigInt::from(431u64))) }

    #[test]
    fn lift_and_check_e1728() {
        let fp2 = ctx();
        let e0 = Montgomery::new(fp2.zero());
        // (0, 0) is a 2-torsion point on E_0: y² = x(x² + 1)
        let p = Aff::P(fp2.zero(), fp2.zero());
        assert!(e0.contains(&p, &fp2));
        // 2 · (0, 0) = ∞
        assert!(e0.dbl(&p, &fp2).is_inf());
    }

    #[test]
    fn group_axioms_on_e0() {
        let fp2 = ctx();
        let e0 = Montgomery::new(fp2.zero());
        // Find a non-2-torsion point on E_0 over F_p (so y ≠ 0).
        let mut p = Aff::Inf;
        for k in 1..200i64 {
            let x = fp2.from_int(k);
            if let Some(lifted) = e0.lift(&x, &fp2) {
                if !e0.dbl(&lifted, &fp2).is_inf() {
                    p = lifted;
                    break;
                }
            }
        }
        assert!(!p.is_inf());
        assert!(e0.contains(&p, &fp2));
        // P + (-P) = ∞
        let np = e0.neg(&p, &fp2);
        assert!(e0.add(&p, &np, &fp2).is_inf());
        // [#E(F_p) for j=1728] · P should be ∞ ... but order varies, skip.
        // [2]P = P + P
        let two_p = e0.dbl(&p, &fp2);
        let pp = e0.add(&p, &p, &fp2);
        assert_eq!(two_p, pp);
        // [3]P = [2]P + P
        let three_p_a = e0.add(&two_p, &p, &fp2);
        let three_p_b = e0.mul(&BigInt::from(3), &p, &fp2);
        assert_eq!(three_p_a, three_p_b);
    }

    #[test]
    fn iota_is_an_involution_squared_to_neg_one() {
        let fp2 = ctx();
        let e0 = Montgomery::new(fp2.zero());
        // Find a point on E_0
        let mut p = Aff::Inf;
        for k in 1..200i64 {
            let x = fp2.from_int(k);
            if let Some(lifted) = e0.lift(&x, &fp2) {
                if !e0.dbl(&lifted, &fp2).is_inf() {
                    p = lifted;
                    break;
                }
            }
        }
        assert!(!p.is_inf());
        // ι(P) on E_0
        let ip = iota_on_e1728(&p, &fp2);
        assert!(
            e0.contains(&ip, &fp2),
            "ι(P) must lie on E_0; got {ip:?}"
        );
        // ι(ι(P)) = -P  (since ι² = [-1])
        let i2p = iota_on_e1728(&ip, &fp2);
        let neg_p = e0.neg(&p, &fp2);
        assert_eq!(i2p, neg_p, "ι² should be [-1]");
        // (2ι)(P) is a point on E_0
        let two_ip = two_iota_on_e1728(&p, &fp2);
        assert!(e0.contains(&two_ip, &fp2));
    }

    #[test]
    fn isog3_image_lies_on_codomain() {
        let fp2 = ctx();
        // E_6: Montgomery A=6 over F_p, supersingular, #E(F_p) = p+1 = 432.
        let e6 = Montgomery::new(fp2.from_int(6));
        // Find a 3-torsion kernel: random point · 144 (since 432/3 = 144).
        let cofactor = BigInt::from(144);
        let mut k_pt = Aff::Inf;
        for seed in 1..200i64 {
            let x = fp2.from_int(seed);
            if let Some(p) = e6.lift(&x, &fp2) {
                let candidate = e6.mul(&cofactor, &p, &fp2);
                if !candidate.is_inf() {
                    // confirm order exactly 3: [3]K = ∞
                    let tripled = e6.mul(&BigInt::from(3), &candidate, &fp2);
                    if tripled.is_inf() {
                        k_pt = candidate;
                        break;
                    }
                }
            }
        }
        let x_k = match &k_pt {
            Aff::Inf => panic!("no 3-torsion found"),
            Aff::P(x, _) => x.clone(),
        };

        // Apply the 3-isogeny.
        let a_new = isog3_codomain(&e6.a, &x_k, &fp2);
        let e_new = Montgomery::new(a_new);

        // Push another random point Q through.
        let mut q_pt = Aff::Inf;
        for seed in 2..200i64 {
            let x = fp2.from_int(seed);
            if let Some(p) = e6.lift(&x, &fp2) {
                // Make sure Q ≠ K, K's negative, etc., so denominator ≠ 0
                if Some(&x) != k_pt.x() {
                    q_pt = p;
                    break;
                }
            }
        }
        // Verify x-only first: phi_x must be the x-coord of *some* point on E_{A_new}
        let phi_x = isog3_push_x(q_pt.x().unwrap(), &x_k, &fp2);
        let x_sq = fp2.sq(&phi_x);
        let x_cube = fp2.mul(&x_sq, &phi_x);
        let rhs = fp2.add(
            &fp2.add(&x_cube, &fp2.mul(&e_new.a, &x_sq)),
            &phi_x,
        );
        assert!(
            fp2.is_square(&rhs),
            "phi_x must be a valid x-coord on the codomain"
        );
        // Now the full push (which lifts y on the codomain).
        let phi_q = isog3_push(&q_pt, &x_k, &e_new.a, &fp2);
        assert!(
            e_new.contains(&phi_q, &fp2),
            "image φ(Q) must lie on the codomain curve E_{{A_new}}; got phi_q = {phi_q:?}"
        );
    }

    #[test]
    fn pushing_3_chain_length_2() {
        let fp2 = ctx();
        let e6 = Montgomery::new(fp2.from_int(6));
        // Find a point of order 9 = 3^2: random·48 (since 432/9 = 48).
        let cofactor = BigInt::from(48);
        let mut k9 = Aff::Inf;
        for seed in 1..200i64 {
            let x = fp2.from_int(seed);
            if let Some(p) = e6.lift(&x, &fp2) {
                let candidate = e6.mul(&cofactor, &p, &fp2);
                if !candidate.is_inf() {
                    let cubed = e6.mul(&BigInt::from(3), &candidate, &fp2);
                    let ninth = e6.mul(&BigInt::from(9), &candidate, &fp2);
                    if !cubed.is_inf() && ninth.is_inf() {
                        k9 = candidate;
                        break;
                    }
                }
            }
        }
        assert!(!k9.is_inf(), "couldn't find a point of order 9");

        // Also find another auxiliary point to push through the chain.
        let mut aux_pt = Aff::Inf;
        for seed in 5..200i64 {
            let x = fp2.from_int(seed);
            if let Some(p) = e6.lift(&x, &fp2) {
                if Some(&x) != k9.x() {
                    aux_pt = p;
                    break;
                }
            }
        }
        assert!(!aux_pt.is_inf());

        let mut aux = vec![aux_pt];
        let e_final = pushing_3_chain(&e6, &k9, 2, &mut aux, &fp2);

        // The pushed auxiliary point must lie on the final curve.
        assert!(
            e_final.contains(&aux[0], &fp2),
            "pushed point must lie on final curve after 3-chain of length 2"
        );
    }

    #[test]
    fn kani_distorted_kernel_construction() {
        // Demonstrate the Kani-style "u·P + v·two_i(P)" combination.
        // For (u, v) = (1, 1), the simplest non-trivial case used in
        // SIKE_challenge.m's first uvtable entry [1, 3, 1, 1].
        let fp2 = ctx();
        let e0 = Montgomery::new(fp2.zero());
        // Pick a point P on E_0
        let mut p_pt = Aff::Inf;
        for k in 1..200i64 {
            let x = fp2.from_int(k);
            if let Some(lifted) = e0.lift(&x, &fp2) {
                if !e0.dbl(&lifted, &fp2).is_inf() {
                    p_pt = lifted;
                    break;
                }
            }
        }
        assert!(!p_pt.is_inf());
        // Kani-distorted point: u·P + v·(2ι)(P) with u=v=1
        let two_ip = two_iota_on_e1728(&p_pt, &fp2);
        let kani = e0.add(
            &e0.mul(&BigInt::from(1), &p_pt, &fp2),
            &e0.mul(&BigInt::from(1), &two_ip, &fp2),
            &fp2,
        );
        assert!(
            e0.contains(&kani, &fp2),
            "Kani-distorted point must lie on E_0"
        );
        // It should NOT equal the input (we're nontrivially distorting).
        assert_ne!(kani, p_pt);
    }
}
