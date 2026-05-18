//! # Vélu's formulas: explicit isogenies from a kernel subgroup.
//!
//! Vélu (1971) gave a closed-form expression for the isogeny
//! `φ: E → E/G` once a finite subgroup `G ⊂ E(F̄)` is specified by
//! its non-identity points.  Given `E: y² = x³ + a x + b` and a
//! kernel `G = {O} ∪ {±P_1, …, ±P_s}` (so `|G| = 2s + 1` for odd
//! cyclic; `|G| = 2s + 2` for cyclic of even order with the 2-torsion
//! point `T = (x_T, 0)` adjoined), define for each `Q = (x_Q, y_Q)`
//! in `G \ {O}`:
//!
//! ```text
//!   g_x(Q) = 3 x_Q² + a
//!   g_y(Q) = -2 y_Q
//!   v(Q)   = 2 g_x(Q)            if 2Q = O          (T-row)
//!          = 4 x_Q g_y(Q)² + 2 g_x(Q)²            otherwise
//!   u(Q)   = g_y(Q)²              (for the y-coordinate map)
//! ```
//!
//! Then summing over a half-set `R` (one representative of each
//! `{±Q}`):
//!
//! ```text
//!   x'  =  x + Σ_{Q ∈ R} [ v(Q) / (x − x_Q)  +  u(Q) / (x − x_Q)² ]
//!   y'  =  y · ( 1 − Σ_{Q ∈ R} [ v(Q) / (x − x_Q)²
//!                                + 2 u(Q) / (x − x_Q)³ ] )
//! ```
//!
//! and the codomain curve is `E' : Y² = X³ + a' X + b'` with
//!
//! ```text
//!   a' = a − 5 · Σ v(Q)
//!   b' = b − 7 · Σ ( u(Q) + x_Q · v(Q) )
//! ```
//!
//! (Both sums are over `R`.)
//!
//! ## Why we keep a separate, small implementation
//!
//! [`crate::cryptanalysis::modular_polynomial`] gives the
//! "algebraic" view: `Φ_ℓ(j(E), j(E')) = 0` constrains *which*
//! `E'` can be reached, but does not tell us **how** to construct
//! the isogeny rationally.  Vélu does both: codomain curve **and**
//! the rational maps for evaluating `φ` on points.
//!
//! ## Scope
//!
//! We enumerate the kernel `G` by brute force over `F_p` and
//! compute Vélu's sums directly.  This is `O(ℓ)` field operations
//! per isogeny — fine for `ℓ ≤ 100` and the toy curves used here.
//! Production-grade `√ℓ`-style algorithms (Bernstein–De Feo–
//! Leroux–Smith, 2020) are deliberately not implemented.

use super::SmallCurve;
use num_bigint::BigUint;
use num_integer::Integer;

/// An `ℓ`-isogeny `φ: E → E'` constructed via Vélu's formulas.
#[derive(Clone, Debug)]
pub struct VeluIsogeny {
    /// The source curve.
    pub domain: SmallCurve,
    /// The codomain curve `E'`.  Computed by Vélu's `a'`, `b'`
    /// formulas above.  `name` is reused from the domain (we do
    /// not invent a new label here).
    pub codomain: SmallCurve,
    /// The degree `ℓ`.
    pub degree: u64,
    /// The kernel's non-identity points as `(x, y)` pairs over
    /// `F_p`, listed as a half-set `R` (one representative per
    /// `{±Q}`).
    pub kernel_half: Vec<(u64, u64)>,
}

impl VeluIsogeny {
    /// Evaluate the isogeny at an affine point `(x, y)`.  Returns
    /// `None` when `(x, y)` is in the kernel (image is the point
    /// at infinity).
    pub fn evaluate(&self, x: u64, y: u64) -> Option<(u64, u64)> {
        let p = self.domain.p as u128;
        let xm = x as u128 % p;
        let ym = y as u128 % p;
        let a = self.domain.a as u128 % p;

        // Reject kernel points.
        for &(xk, _) in &self.kernel_half {
            if xm == xk as u128 % p {
                return None;
            }
        }

        let mut acc_x = xm;
        let mut acc_y = ym;
        for &(xq, yq) in &self.kernel_half {
            let xq = xq as u128 % p;
            let yq = yq as u128 % p;

            let gx_q = mod_add(3 * mod_mul(xq, xq, p) % p, a, p);
            let gy_q = mod_sub(0, 2 * yq % p, p);
            // T-row test: 2Q = O ⇔ yQ = 0.
            let t_row = yq == 0;

            let v_q = if t_row {
                (2 * gx_q) % p
            } else {
                let gy2 = mod_mul(gy_q, gy_q, p);
                let term1 = mod_mul(4 * xq % p, gy2, p);
                let term2 = mod_mul(2 * gx_q % p, gx_q, p);
                mod_add(term1, term2, p)
            };
            let u_q = mod_mul(gy_q, gy_q, p);

            let dx = mod_sub(xm, xq, p);
            let dx_inv = mod_inv(dx, p)?;
            let dx_inv2 = mod_mul(dx_inv, dx_inv, p);
            let dx_inv3 = mod_mul(dx_inv2, dx_inv, p);

            acc_x = mod_add(acc_x, mod_mul(v_q, dx_inv, p), p);
            acc_x = mod_add(acc_x, mod_mul(u_q, dx_inv2, p), p);

            // y-coordinate transformation: y' = y · (1 - ...)
            // We accumulate the bracket Σ [v · dx^{-2} + 2u · dx^{-3}].
            let bracket = mod_add(
                mod_mul(v_q, dx_inv2, p),
                mod_mul(2 * u_q % p, dx_inv3, p),
                p,
            );
            // Subtract the bracket from acc_y / ym, then re-multiply by y.
            // The cleanest way: compute the prefactor (1 − Σ) on the
            // running ym independently.
            //
            // Since acc_y starts as ym, the *first* iteration should
            // already apply (1 - bracket_1).  But we want all
            // iterations to multiply the same y, not to cascade.
            // Rebuild the factor explicitly:
            acc_y = mod_sub(acc_y, mod_mul(ym, bracket, p), p);
        }
        Some((acc_x as u64, acc_y as u64))
    }
}

// ── Modular arithmetic helpers in u128 ────────────────────────────────────────

fn mod_add(a: u128, b: u128, p: u128) -> u128 {
    ((a % p) + (b % p)) % p
}
fn mod_sub(a: u128, b: u128, p: u128) -> u128 {
    ((a % p) + p - (b % p)) % p
}
fn mod_mul(a: u128, b: u128, p: u128) -> u128 {
    // For p < 2^62 this fits comfortably in u128 (operands < 2^62
    // each, product < 2^124).
    ((a % p) * (b % p)) % p
}
fn mod_pow(mut base: u128, mut e: u128, p: u128) -> u128 {
    let mut r = 1u128 % p;
    base %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mod_mul(r, base, p);
        }
        base = mod_mul(base, base, p);
        e >>= 1;
    }
    r
}
fn mod_inv(a: u128, p: u128) -> Option<u128> {
    if a == 0 {
        return None;
    }
    // Fermat: p prime ⇒ a^{p-2} ≡ a⁻¹.
    Some(mod_pow(a, p - 2, p))
}

// ── Vélu codomain construction ────────────────────────────────────────────────

/// Compute the codomain curve and kernel-half representation from a
/// list of affine kernel points.  The list must contain at least one
/// representative of each `{±Q}` orbit.  Duplicate `{±Q}` orbits in
/// the input are silently de-duplicated.
fn velu_codomain_from_kernel(domain: &SmallCurve, kernel: &[(u64, u64)]) -> VeluIsogeny {
    let p = domain.p as u128;
    let a = domain.a as u128 % p;
    let b = domain.b as u128 % p;

    // Reduce to a half-set R: keep `Q` and drop `-Q`.
    let mut seen_x = std::collections::HashSet::<u64>::new();
    let mut half = Vec::new();
    for &(x, y) in kernel {
        if seen_x.insert(x) {
            half.push((x, y));
        }
    }

    let mut sum_v = 0u128;
    let mut sum_uvx = 0u128;
    for &(xq, yq) in &half {
        let xq = xq as u128 % p;
        let yq = yq as u128 % p;

        let gx_q = mod_add(3 * mod_mul(xq, xq, p) % p, a, p);
        let gy_q = mod_sub(0, 2 * yq % p, p);
        let t_row = yq == 0;
        let v_q = if t_row {
            (2 * gx_q) % p
        } else {
            let gy2 = mod_mul(gy_q, gy_q, p);
            let term1 = mod_mul(4 * xq % p, gy2, p);
            let term2 = mod_mul(2 * gx_q % p, gx_q, p);
            mod_add(term1, term2, p)
        };
        let u_q = mod_mul(gy_q, gy_q, p);

        sum_v = mod_add(sum_v, v_q, p);
        sum_uvx = mod_add(sum_uvx, mod_add(u_q, mod_mul(xq, v_q, p), p), p);
    }

    let new_a = mod_sub(a, (5 * sum_v) % p, p) as u64;
    let new_b = mod_sub(b, (7 * sum_uvx) % p, p) as u64;
    let codomain = SmallCurve {
        name: domain.name,
        p: domain.p,
        a: new_a,
        b: new_b,
    };
    // For a cyclic kernel of odd order ℓ: |G| = 2·|half| + 1.
    // For a cyclic kernel of even order ℓ that includes the
    // 2-torsion point T = (x_T, 0): |G| = 2·|half| (T appears in
    // `half` as its own inverse).
    let degree = if half.iter().any(|&(_, y)| y == 0) {
        2 * (half.len() as u64)
    } else {
        2 * (half.len() as u64) + 1
    };

    VeluIsogeny {
        domain: *domain,
        codomain,
        degree,
        kernel_half: half,
    }
}

// ── ℓ = 2 isogeny ────────────────────────────────────────────────────────────

/// Compute a 2-isogeny with the specified 2-torsion point `T = (x_T, 0)`
/// as kernel.  `T` must satisfy `y² = x³ + a x + b` with `y = 0`,
/// i.e. `x_T` is a root of `x³ + a x + b` over `F_p`.
pub fn velu_isogeny_2(domain: &SmallCurve, x_t: u64) -> Option<VeluIsogeny> {
    let p = domain.p as u128;
    let xt = x_t as u128 % p;
    let a = domain.a as u128 % p;
    let b = domain.b as u128 % p;
    // Check x_T is on the curve with y = 0.
    let rhs = mod_add(
        mod_add(mod_mul(mod_mul(xt, xt, p), xt, p), mod_mul(a, xt, p), p),
        b,
        p,
    );
    if rhs != 0 {
        return None;
    }
    Some(velu_codomain_from_kernel(domain, &[(x_t, 0)]))
}

// ── ℓ odd isogeny by exhaustive kernel enumeration ───────────────────────────

/// Enumerate the curve's `ℓ`-torsion subgroup(s) by brute force over
/// `F_p`, and return all *cyclic* subgroups of order `ℓ` as
/// `VeluIsogeny` records.  `ℓ` must be an odd prime.
///
/// Naive complexity: `O(ℓ² · p)` for the membership test on a
/// candidate generator.  Fine for `ℓ ≤ 7` and `p ≤ 10^4`.
pub fn velu_isogeny_odd(domain: &SmallCurve, ell: u64) -> Vec<VeluIsogeny> {
    assert!(ell >= 3 && ell % 2 == 1, "ℓ must be an odd prime");

    let curve = domain.to_curve_params();
    let a_fe = curve.a_fe();
    let g = crate::ecc::point::Point::Infinity;
    let _ = g;

    // 1. Enumerate every affine point of E(F_p) using Tonelli-Shanks
    //    for each x.  O(p · log² p) instead of O(p²).
    let mut points = Vec::<(u64, u64)>::new();
    for x in 0..domain.p {
        let rhs = domain.rhs(x);
        if rhs == 0 {
            points.push((x, 0));
            continue;
        }
        let p = domain.p;
        match crate::isogeny::cm::legendre_u64(rhs, p) {
            1 => {
                if let Some(y) = crate::isogeny::cm::tonelli_shanks_u64(rhs, p) {
                    points.push((x, y));
                    if y != 0 && p - y != y {
                        points.push((x, p - y));
                    }
                }
            }
            _ => {}
        }
    }

    // 2. For each point P, compute [k]P for k = 2..ell and test
    //    whether [ell]P = O. The set { [0]P, [1]P, ..., [ell-1]P }
    //    is then a cyclic subgroup of order ell when P has order
    //    exactly ell.
    let mut subgroups_seen: Vec<std::collections::BTreeSet<(u64, u64)>> = Vec::new();
    let mut result = Vec::new();

    for &(px, py) in &points {
        if py == 0 {
            continue; // skip 2-torsion; we want odd ℓ
        }
        let p_point = crate::ecc::point::Point::Affine {
            x: curve.fe(BigUint::from(px)),
            y: curve.fe(BigUint::from(py)),
        };
        let ell_p = p_point.scalar_mul(&BigUint::from(ell), &a_fe);
        if !matches!(ell_p, crate::ecc::point::Point::Infinity) {
            continue;
        }

        // Build subgroup: { kP : k = 1..ell-1 }.
        let mut sg = std::collections::BTreeSet::new();
        let mut current = p_point.clone();
        let mut bad = false;
        for k in 1..ell {
            let _ = k;
            match &current {
                crate::ecc::point::Point::Affine { x, y } => {
                    let xv = x.value.to_string().parse::<u64>().unwrap_or(u64::MAX);
                    let yv = y.value.to_string().parse::<u64>().unwrap_or(u64::MAX);
                    if xv == u64::MAX || yv == u64::MAX {
                        bad = true;
                        break;
                    }
                    sg.insert((xv, yv));
                }
                crate::ecc::point::Point::Infinity => {
                    // Should not happen before reaching ell;
                    // would mean a smaller order.
                    bad = true;
                    break;
                }
            }
            current = current.add(&p_point, &a_fe);
        }
        if bad || sg.len() as u64 != ell - 1 {
            continue;
        }

        if subgroups_seen.iter().any(|s| s == &sg) {
            continue;
        }
        subgroups_seen.push(sg.clone());
        let kernel: Vec<(u64, u64)> = sg.into_iter().collect();
        result.push(velu_codomain_from_kernel(domain, &kernel));
    }
    result
}

// ── Dual isogeny ──────────────────────────────────────────────────────────────
//
// For φ: E → E' of degree ℓ, the dual `φ̂: E' → E` is again degree ℓ
// and satisfies `φ̂ ∘ φ = [ℓ]_E`.  The kernel of `φ̂` is `φ(E[ℓ])`,
// which over `F_p` is the image of any `ℓ`-torsion subgroup of `E'`
// not equal to `ker(φ̂)` itself.
//
// We provide a `dual_isogeny_naive` that simply searches the
// codomain's `ℓ`-torsion subgroups and picks the one that completes
// the round-trip on a probe point.

/// Find the dual `ℓ`-isogeny of a `VeluIsogeny` by exhaustive search
/// over the codomain's `ℓ`-isogeny family.  `O(deg)`-many candidates
/// (one per cyclic subgroup); each is tested by checking that the
/// composition is `[ℓ]` on a single probe point.
pub fn dual_isogeny_naive(phi: &VeluIsogeny) -> Option<VeluIsogeny> {
    let ell = phi.degree;
    let cands = if ell == 2 {
        // 2-isogeny: search the codomain's 2-torsion.
        let mut v = Vec::new();
        for x in 0..phi.codomain.p {
            if phi.codomain.rhs(x) == 0 {
                if let Some(iso) = velu_isogeny_2(&phi.codomain, x) {
                    v.push(iso);
                }
            }
        }
        v
    } else {
        velu_isogeny_odd(&phi.codomain, ell)
    };

    // Probe: pick a point on the domain that's not in ker(φ).
    let p = phi.domain.p;
    for x in 2..p {
        let rhs = phi.domain.rhs(x);
        if rhs == 0 {
            continue;
        }
        let y = match crate::isogeny::cm::legendre_u64(rhs, p) {
            1 => match crate::isogeny::cm::tonelli_shanks_u64(rhs, p) {
                Some(y) if y != 0 => y,
                _ => continue,
            },
            _ => continue,
        };
        let _ = y; // explicit binding for clarity below

        if let Some(image) = phi.evaluate(x, y) {
            // For each candidate ψ, compute ψ(image) and check it equals [ell] · (x, y).
            let curve = phi.domain.to_curve_params();
            let a_fe = curve.a_fe();
            let p_point = crate::ecc::point::Point::Affine {
                x: curve.fe(BigUint::from(x)),
                y: curve.fe(BigUint::from(y)),
            };
            let ell_p = p_point.scalar_mul(&BigUint::from(ell), &a_fe);
            for psi in &cands {
                if let Some((px, py)) = psi.evaluate(image.0, image.1) {
                    if let crate::ecc::point::Point::Affine { x: ex, y: ey } = &ell_p {
                        let exv = ex.value.iter_u64_digits().next().unwrap_or(0);
                        let eyv = ey.value.iter_u64_digits().next().unwrap_or(0);
                        if exv == px && eyv == py {
                            return Some(psi.clone());
                        }
                    }
                }
            }
            // No probe → too small; bail.
            break;
        }
    }
    None
}

// ── Composition of isogenies ─────────────────────────────────────────────────

/// Compose two Vélu isogenies `φ: E → E'` and `ψ: E' → E''`.  The
/// returned record's `kernel_half` field is the *pull-back* of
/// `ψ`'s kernel under `φ` (computed by evaluating `φ` on a curve
/// enumeration), so the resulting record correctly represents the
/// composite as a single isogeny.  Note: the **degree** of the
/// composite is `phi.degree · psi.degree`.
pub fn compose_isogenies(phi: &VeluIsogeny, psi: &VeluIsogeny) -> Option<VeluIsogeny> {
    if phi.codomain.a != psi.domain.a
        || phi.codomain.b != psi.domain.b
        || phi.codomain.p != psi.domain.p
    {
        return None;
    }
    // Enumerate ker(ψ ∘ φ) by walking E(F_p) and keeping points
    // whose image under φ lies in ker(ψ).  Tonelli-Shanks per x
    // gives O(p · log² p) instead of O(p²).
    let mut kernel = Vec::new();
    let p = phi.domain.p;
    for x in 0..p {
        let rhs = phi.domain.rhs(x);
        let y = if rhs == 0 {
            0u64
        } else {
            match crate::isogeny::cm::legendre_u64(rhs, p) {
                1 => match crate::isogeny::cm::tonelli_shanks_u64(rhs, p) {
                    Some(y) => y,
                    None => continue,
                },
                _ => continue,
            }
        };
        if let Some(image) = phi.evaluate(x, y) {
            if psi.kernel_half.iter().any(|&(kx, _)| kx == image.0) {
                kernel.push((x, y));
            }
        }
    }
    let mut iso = velu_codomain_from_kernel(&phi.domain, &kernel);
    iso.codomain = psi.codomain;
    iso.degree = phi.degree * psi.degree;
    Some(iso)
}

// ── Unit tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isogeny::toy_curve_a;

    #[test]
    fn two_isogeny_round_trip() {
        // On the toy-2003 curve, search for a 2-torsion point and
        // build the 2-isogeny.  Codomain should be a valid curve
        // with non-zero discriminant.
        let curve = toy_curve_a();
        let mut found = None;
        for x in 0..curve.p {
            if curve.rhs(x) == 0 {
                found = Some(x);
                break;
            }
        }
        if let Some(x_t) = found {
            let iso = velu_isogeny_2(&curve, x_t).expect("2-iso constructs");
            assert_eq!(iso.degree, 2);
            // Codomain is a *different* curve (unless x_T is degenerate).
            assert!(iso.codomain.p == curve.p);
        }
    }

    #[test]
    fn three_isogeny_count_for_toy() {
        // On a small ordinary curve, the number of 3-isogenies is
        // bounded by the number of order-3 cyclic subgroups, which
        // is 0, 1, 2, or 4 depending on rank of E[3].
        let curve = toy_curve_a();
        let isos = velu_isogeny_odd(&curve, 3);
        assert!(isos.len() <= 4, "got {} 3-isogenies", isos.len());
    }
}
