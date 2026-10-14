//! `O(ℓ)` Vélu for odd prime `ℓ` on x-only Montgomery curves: images by the
//! Costello–Hisil formula, codomain by the Meyer–Reith twisted-Edwards
//! formula.
//!
//! For a kernel `⟨P⟩` of odd prime order `ℓ`, write `d = (ℓ−1)/2` and
//! `(X_i : Z_i) = x([i]P)`, `1 ≤ i ≤ d`.  Then
//!
//! * `x(φ(Q)) = x · Π_i ((x X_i − Z_i)/(X_i − x Z_i))²`;
//! * with `a_E = A + 2`, `d_E = A − 2`:
//!   `a' = a_E^ℓ · Π(X_i + Z_i)^8`, `d' = d_E^ℓ · Π(X_i − Z_i)^8`, and
//!   `A' = 2(a' + d')/(a' − d')`.
//!
//! Cost: `≈ 6d` `M` for the multiples, `2d + O(log ℓ)` `M` plus one `I`
//! for the codomain, `4d + 4` `M` plus one `I` per image.

use super::counted_fp::{OpCounter, PrimeField};
use super::montgomery::{xadd, xdbl, MontgomeryCurve, XZ};
use super::IsogenyEvalError;

/// The result of evaluating an isogeny: its codomain and the images of the
/// requested points.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyEvaluation {
    /// The codomain's Montgomery coefficient `A'`.
    pub codomain: MontgomeryCurve,
    /// `x(φ(Q))` for each input `x(Q)`, in order; `None` when `Q` lies in
    /// the kernel (its image is the point at infinity).
    pub images: Vec<Option<u64>>,
}

/// Check `ell` is an odd prime.
pub(crate) fn check_degree(ell: u64) -> Result<(), IsogenyEvalError> {
    if ell < 3 || !super::counted_fp::is_prime_u64(ell) {
        return Err(IsogenyEvalError::InvalidDegree(ell));
    }
    Ok(())
}

/// The Edwards-form codomain `A' = 2(a' + d')/(a' − d')` from
/// `prod_plus = c·Π(x_s + 1)` and `prod_minus = c'·Π(x_s − 1)` with
/// `c^8 = c'^8` (eighth powers cancel common factors and signs).
pub(crate) fn edwards_codomain<C: OpCounter>(
    f: &PrimeField<C>,
    curve: &MontgomeryCurve,
    ell: u64,
    prod_plus: u64,
    prod_minus: u64,
) -> Result<MontgomeryCurve, IsogenyEvalError> {
    let pow8 = |v: u64| {
        let v2 = f.sqr(v);
        let v4 = f.sqr(v2);
        f.sqr(v4)
    };
    let a = curve.a;
    let ae = f.mul(f.pow(f.add(a, 2), ell as u128), pow8(prod_plus));
    let de = f.mul(f.pow(f.sub(a, 2), ell as u128), pow8(prod_minus));
    let num = f.mul_small(2, f.add(ae, de));
    let a_new = f
        .div(num, f.sub(ae, de))
        .ok_or(IsogenyEvalError::SingularCurve)?;
    MontgomeryCurve::new(f, a_new)
}

/// `x([i]P)` for `1 ≤ i ≤ d` (`d ≥ 1`), by one doubling and `d − 2`
/// differential additions.
pub(crate) fn first_multiples<C: OpCounter>(
    f: &PrimeField<C>,
    curve: &MontgomeryCurve,
    xp: u64,
    d: usize,
) -> Vec<XZ<u64>> {
    let p1 = XZ::affine(xp);
    let mut mult = Vec::with_capacity(d);
    mult.push(p1.clone());
    if d >= 2 {
        mult.push(xdbl(f, &p1, &curve.a24));
    }
    for i in 2..d {
        let next = xadd(f, &mult[i - 1], &p1, &mult[i - 2]);
        mult.push(next);
    }
    mult
}

/// Evaluate the `ℓ`-isogeny with kernel `⟨P⟩`, `x(P) = kernel_x`, at the
/// x-coordinates `points` (`O(ℓ)` per point).
///
/// Preconditions: `ell` is an odd prime and `kernel_x` is the
/// x-coordinate of a point of exact order `ell` on the curve or its twist
/// (check with [`super::montgomery::is_kernel_generator`]); the points are
/// canonical field elements.  A wrong kernel gives a meaningless (but
/// non-panicking) result.
///
/// Errors: [`IsogenyEvalError::InvalidDegree`] if `ell` is not an odd
/// prime, [`IsogenyEvalError::SingularCurve`] if the computed codomain is
/// singular (impossible for a valid kernel).
pub fn velu<C: OpCounter>(
    f: &PrimeField<C>,
    curve: &MontgomeryCurve,
    ell: u64,
    kernel_x: u64,
    points: &[u64],
) -> Result<IsogenyEvaluation, IsogenyEvalError> {
    check_degree(ell)?;
    let d = ((ell - 1) / 2) as usize;
    let mult = first_multiples(f, curve, kernel_x, d);
    // (X_i − Z_i, X_i + Z_i).
    let sums: Vec<(u64, u64)> = mult
        .iter()
        .map(|q| (f.sub(q.x, q.z), f.add(q.x, q.z)))
        .collect();
    let (mut pp, mut pm) = (1u64, 1u64);
    for &(m, p) in &sums {
        pm = f.mul(pm, m);
        pp = f.mul(pp, p);
    }
    let codomain = edwards_codomain(f, curve, ell, pp, pm)?;
    let images = points
        .iter()
        .map(|&x| {
            let (xm, xp) = (f.sub(x, 1), f.add(x, 1));
            let (mut num, mut den) = (1u64, 1u64);
            for &(m, p) in &sums {
                let t0 = f.mul(m, xp);
                let t1 = f.mul(p, xm);
                num = f.mul(num, f.add(t0, t1));
                den = f.mul(den, f.sub(t0, t1));
            }
            let xn = f.mul(x, f.sqr(num));
            let zn = f.sqr(den);
            XZ { x: xn, z: zn }.to_affine(f)
        })
        .collect();
    Ok(IsogenyEvaluation { codomain, images })
}

/// The kernel polynomial `Π_{i=1}^{d} (Z − x([i]P))` on the Montgomery
/// x-line (monic, degree `(ℓ−1)/2`).  Same preconditions as [`velu`];
/// errors with [`IsogenyEvalError::NotAKernelPoint`] if some multiple
/// `[i]P`, `i ≤ d`, is the point at infinity.
pub fn kernel_polynomial<C: OpCounter>(
    f: &PrimeField<C>,
    curve: &MontgomeryCurve,
    ell: u64,
    kernel_x: u64,
) -> Result<Vec<u64>, IsogenyEvalError> {
    check_degree(ell)?;
    let d = ((ell - 1) / 2) as usize;
    let mult = first_multiples(f, curve, kernel_x, d);
    let xs =
        super::montgomery::batch_to_affine(f, &mult).ok_or(IsogenyEvalError::NotAKernelPoint)?;
    let leaves = xs.into_iter().map(|x| vec![f.neg(x), 1]).collect();
    Ok(super::poly::product(f, leaves))
}
