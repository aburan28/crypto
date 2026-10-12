//! Vélu isogenies of degree 2 and 3 on short Weierstrass curves over `F_p`,
//! with full `(x, y)` point evaluation, `F_p`-isomorphisms, and a certified
//! dual.
//!
//! The library's [`crate::cryptanalysis::isogeny_walk::kernel`] gives
//! codomains of odd degree only, without point evaluation or duals.  This
//! module covers `ℓ ∈ {2, 3}` completely:
//!
//! * kernels: the rational roots `x₀` of `x³ + ax + b` (`ℓ = 2`) or of
//!   `ψ₃ = 3x⁴ + 6ax² + 12bx − a²` (`ℓ = 3`; the kernel `{O, ±Q}` is
//!   `F_p`-rational even when `y(Q)` is not);
//! * Vélu: `v = 3x₀² + a`, `u = 0` (`ℓ = 2`) or `v = 2(3x₀² + a)`,
//!   `u = 4(x₀³ + ax₀ + b)` (`ℓ = 3`); `w = u + x₀v`; codomain
//!   `a' = a − 5v`, `b' = b − 7w`; with `z = 1/(x − x₀)`,
//!   `φ(x, y) = (x + vz + uz², y(1 − vz² − 2uz³))`;
//! * the dual: a Vélu isogeny on the codomain whose codomain has the
//!   domain's `j`, followed by an `F_p`-isomorphism `(x, y) ↦ (u²x, u³y)`,
//!   accepted only if `ψ(φ(Q)) = [ℓ]Q` for a point `Q` of order `> 4ℓ²`.
//!   That check is a proof: `ψφ − [ℓ]` has degree `≤ (ℓ + ℓ)² = 4ℓ²`, so if
//!   it is nonzero its kernel cannot contain `⟨Q⟩`.
//!
//! Point evaluation costs `1 I + 7 M` (plus `2 M` for the isomorphism).
//! The kernel search ([`kernel_x_coordinates`]) and the isomorphism search
//! ([`isomorphisms`]) use the library's root finder over its own field and
//! are not counted.

use num_bigint::BigUint;
use num_traits::ToPrimitive;

use super::counted_fp::{OpCounter, PrimeField};
use super::IsogenyEvalError;
use crate::cryptanalysis::isogeny_walk::field::Field as LibField;
use crate::cryptanalysis::isogeny_walk::poly as lib_poly;

/// An affine point, `None` for the point at infinity.
pub type WPoint = Option<(u64, u64)>;

/// A short Weierstrass curve `y² = x³ + ax + b` over `F_p` (`p ≥ 5`),
/// nonsingular.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShortWeierstrass {
    /// `a`.
    pub a: u64,
    /// `b`.
    pub b: u64,
}

impl ShortWeierstrass {
    /// The curve `y² = x³ + ax + b` (coefficients reduced mod `p`); errors
    /// with [`IsogenyEvalError::SingularCurve`] if `4a³ + 27b² = 0`.
    pub fn new<C: OpCounter>(f: &PrimeField<C>, a: u64, b: u64) -> Result<Self, IsogenyEvalError> {
        let e = ShortWeierstrass {
            a: f.elem(a),
            b: f.elem(b),
        };
        if e.discriminant_part(f) == 0 {
            return Err(IsogenyEvalError::SingularCurve);
        }
        Ok(e)
    }

    /// `4a³ + 27b²`.
    fn discriminant_part<C: OpCounter>(&self, f: &PrimeField<C>) -> u64 {
        let a3 = f.mul(f.sqr(self.a), self.a);
        f.add(f.mul_small(4, a3), f.mul_small(27, f.sqr(self.b)))
    }

    /// `x³ + ax + b`.
    pub fn rhs<C: OpCounter>(&self, f: &PrimeField<C>, x: u64) -> u64 {
        f.add(f.mul(f.add(f.sqr(x), self.a), x), self.b)
    }

    /// Whether `p` lies on the curve (the point at infinity does).
    pub fn contains<C: OpCounter>(&self, f: &PrimeField<C>, p: WPoint) -> bool {
        match p {
            None => true,
            Some((x, y)) => x < f.p() && y < f.p() && f.sqr(y) == self.rhs(f, x),
        }
    }

    /// `j = 1728 · 4a³ / (4a³ + 27b²)`.
    pub fn j_invariant<C: OpCounter>(&self, f: &PrimeField<C>) -> u64 {
        let a3 = f.mul(f.sqr(self.a), self.a);
        let num = f.mul_small(1728 * 4, a3);
        // Nonsingular by construction, so the denominator is invertible.
        f.div(num, self.discriminant_part(f)).unwrap_or(0)
    }

    /// `−P`.
    pub fn neg<C: OpCounter>(&self, f: &PrimeField<C>, p: WPoint) -> WPoint {
        p.map(|(x, y)| (x, f.neg(y)))
    }

    /// `P + Q` (affine; one `I` and two or three `M` plus the squaring).
    pub fn add<C: OpCounter>(&self, f: &PrimeField<C>, p: WPoint, q: WPoint) -> WPoint {
        let (x1, y1) = match p {
            None => return q,
            Some(v) => v,
        };
        let (x2, y2) = match q {
            None => return p,
            Some(v) => v,
        };
        let lam = if x1 == x2 {
            if y1 != y2 || y1 == 0 {
                return None;
            }
            let num = f.add(f.mul_small(3, f.sqr(x1)), self.a);
            f.div(num, f.add(y1, y1))?
        } else {
            f.div(f.sub(y2, y1), f.sub(x2, x1))?
        };
        let x3 = f.sub(f.sub(f.sqr(lam), x1), x2);
        let y3 = f.sub(f.mul(lam, f.sub(x1, x3)), y1);
        Some((x3, y3))
    }

    /// `[k]P` by double-and-add.
    pub fn mul<C: OpCounter>(&self, f: &PrimeField<C>, p: WPoint, k: u128) -> WPoint {
        let mut acc: WPoint = None;
        let mut base = p;
        let mut k = k;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(f, acc, base);
            }
            k >>= 1;
            if k > 0 {
                base = self.add(f, base, base);
            }
        }
        acc
    }

    /// The point with the least `x ≥ start` (and the `y` returned by
    /// [`PrimeField::sqrt`]), scanning at most `limit` values of `x`.
    pub fn point_from_x<C: OpCounter>(
        &self,
        f: &PrimeField<C>,
        start: u64,
        limit: u64,
    ) -> Option<(u64, u64)> {
        (0..limit)
            .map(|i| f.elem(start.wrapping_add(i)))
            .find_map(|x| f.sqrt(self.rhs(f, x)).map(|y| (x, y)))
    }
}

fn lib_field(p: u64) -> Option<LibField> {
    LibField::new(&BigUint::from(p))
}

fn lib_roots(f: &LibField, coeffs: &[u64]) -> Vec<u64> {
    let poly: Vec<_> = coeffs.iter().map(|&c| f.from_u64(c)).collect();
    lib_poly::roots(f, &poly, 7)
        .iter()
        .filter_map(|r| f.to_big(r).to_u64())
        .collect()
}

/// The `x`-coordinates `x₀` of the `F_p`-rational kernels of degree
/// `ell ∈ {2, 3}` on `e`, sorted.  Errors with
/// [`IsogenyEvalError::InvalidDegree`] for other `ell`.
pub fn kernel_x_coordinates<C: OpCounter>(
    f: &PrimeField<C>,
    e: &ShortWeierstrass,
    ell: u64,
) -> Result<Vec<u64>, IsogenyEvalError> {
    let coeffs = match ell {
        2 => vec![e.b, e.a, 0, 1],
        3 => vec![
            f.neg(f.sqr(e.a)),
            f.mul_small(12, e.b),
            f.mul_small(6, e.a),
            0,
            3,
        ],
        _ => return Err(IsogenyEvalError::InvalidDegree(ell)),
    };
    let lf = lib_field(f.p()).ok_or(IsogenyEvalError::SingularCurve)?;
    Ok(lib_roots(&lf, &coeffs))
}

/// Every `u ∈ F_p^*` with `to.a = u⁴·from.a` and `to.b = u⁶·from.b`, i.e.
/// every `F_p`-isomorphism `(x, y) ↦ (u²x, u³y)` from `from` to `to`
/// (two for `j ∉ {0, 1728}`, up to four or six otherwise).  Sorted.
pub fn isomorphisms<C: OpCounter>(
    f: &PrimeField<C>,
    from: &ShortWeierstrass,
    to: &ShortWeierstrass,
) -> Vec<u64> {
    if (from.a == 0) != (to.a == 0) || (from.b == 0) != (to.b == 0) {
        return Vec::new();
    }
    let Some(lf) = lib_field(f.p()) else {
        return Vec::new();
    };
    let candidates = if from.a != 0 {
        let r = f.div(to.a, from.a).unwrap_or(0);
        lib_roots(&lf, &[f.neg(r), 0, 0, 0, 1])
    } else {
        let r = f.div(to.b, from.b).unwrap_or(0);
        lib_roots(&lf, &[f.neg(r), 0, 0, 0, 0, 0, 1])
    };
    candidates
        .into_iter()
        .filter(|&u| {
            let u2 = f.sqr(u);
            let u4 = f.sqr(u2);
            let u6 = f.mul(u4, u2);
            u != 0 && f.mul(u4, from.a) == to.a && f.mul(u6, from.b) == to.b
        })
        .collect()
}

/// A Vélu isogeny of degree 2 or 3 with an `F_p`-rational kernel, optionally
/// followed by an isomorphism `(X, Y) ↦ (u²X, u³Y)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SmallIsogeny {
    ell: u64,
    kernel_x: u64,
    v: u64,
    u: u64,
    domain: ShortWeierstrass,
    codomain: ShortWeierstrass,
    /// `(u², u³)` of the trailing isomorphism, `None` for the identity.
    iso: Option<(u64, u64)>,
}

impl SmallIsogeny {
    /// Vélu's isogeny from `domain` with kernel x-coordinate `kernel_x`.
    ///
    /// Errors: [`IsogenyEvalError::InvalidDegree`] unless `ell ∈ {2, 3}`;
    /// [`IsogenyEvalError::NotAKernelPoint`] unless `kernel_x` is a root of
    /// `x³ + ax + b` (`ℓ = 2`) or `ψ₃` (`ℓ = 3`);
    /// [`IsogenyEvalError::SingularCurve`] if the codomain is singular.
    pub fn new<C: OpCounter>(
        f: &PrimeField<C>,
        domain: &ShortWeierstrass,
        ell: u64,
        kernel_x: u64,
    ) -> Result<Self, IsogenyEvalError> {
        let x0 = f.elem(kernel_x);
        let gx = f.add(f.mul_small(3, f.sqr(x0)), domain.a);
        let (v, u) = match ell {
            2 => {
                if domain.rhs(f, x0) != 0 {
                    return Err(IsogenyEvalError::NotAKernelPoint);
                }
                (gx, 0)
            }
            3 => {
                let x2 = f.sqr(x0);
                let psi3 = f.sub(
                    f.add(
                        f.add(
                            f.mul_small(3, f.sqr(x2)),
                            f.mul_small(6, f.mul(domain.a, x2)),
                        ),
                        f.mul_small(12, f.mul(domain.b, x0)),
                    ),
                    f.sqr(domain.a),
                );
                if psi3 != 0 {
                    return Err(IsogenyEvalError::NotAKernelPoint);
                }
                (f.add(gx, gx), f.mul_small(4, domain.rhs(f, x0)))
            }
            _ => return Err(IsogenyEvalError::InvalidDegree(ell)),
        };
        let w = f.add(u, f.mul(x0, v));
        let codomain = ShortWeierstrass::new(
            f,
            f.sub(domain.a, f.mul_small(5, v)),
            f.sub(domain.b, f.mul_small(7, w)),
        )?;
        Ok(SmallIsogeny {
            ell,
            kernel_x: x0,
            v,
            u,
            domain: *domain,
            codomain,
            iso: None,
        })
    }

    /// The degree `ℓ`.
    pub fn degree(&self) -> u64 {
        self.ell
    }

    /// The kernel x-coordinate `x₀`.
    pub fn kernel_x(&self) -> u64 {
        self.kernel_x
    }

    /// The domain.
    pub fn domain(&self) -> &ShortWeierstrass {
        &self.domain
    }

    /// The codomain (after the isomorphism, if any).
    pub fn codomain(&self) -> &ShortWeierstrass {
        &self.codomain
    }

    /// Vélu's raw codomain coefficients before any isomorphism:
    /// `(a − 5v, b − 7w)`.
    pub fn velu_codomain<C: OpCounter>(&self, f: &PrimeField<C>) -> ShortWeierstrass {
        let w = f.add(self.u, f.mul(self.kernel_x, self.v));
        ShortWeierstrass {
            a: f.sub(self.domain.a, f.mul_small(5, self.v)),
            b: f.sub(self.domain.b, f.mul_small(7, w)),
        }
    }

    /// `(u², u³)` of the trailing isomorphism, `None` if there is none.
    pub fn isomorphism(&self) -> Option<(u64, u64)> {
        self.iso
    }

    /// `φ(P)`.  Kernel points and `O` map to `O`.  `P` must lie on the
    /// domain (not checked).
    pub fn eval<C: OpCounter>(&self, f: &PrimeField<C>, p: WPoint) -> WPoint {
        let (x, y) = p?;
        let z = f.inv(f.sub(x, self.kernel_x))?;
        let z2 = f.sqr(z);
        let z3 = f.mul(z2, z);
        let xx = f.add(f.add(x, f.mul(self.v, z)), f.mul(self.u, z2));
        let d = f.sub(
            f.sub(1, f.mul(self.v, z2)),
            f.mul(f.add(self.u, self.u), z3),
        );
        let yy = f.mul(y, d);
        match self.iso {
            None => Some((xx, yy)),
            Some((s, t)) => Some((f.mul(s, xx), f.mul(t, yy))),
        }
    }

    /// A point of the domain of order `> 4ℓ²`, by scanning `x = 0, 1, …`
    /// (at most `2^16` values).
    fn certification_point<C: OpCounter>(&self, f: &PrimeField<C>) -> Option<(u64, u64)> {
        let e = &self.domain;
        let bound = 4 * self.ell * self.ell;
        let mut x = 0u64;
        while x < (1 << 16) && x < f.p() {
            let q = e.point_from_x(f, x, (1 << 16) - x)?;
            let mut acc = Some(q);
            let mut big = true;
            for _ in 1..bound {
                acc = e.add(f, acc, Some(q));
                if acc.is_none() {
                    big = false;
                    break;
                }
            }
            if big {
                return Some(q);
            }
            x = q.0 + 1;
        }
        None
    }

    /// The dual `φ̂ : codomain → domain`, with `φ̂ ∘ φ = [ℓ]` exactly (same
    /// domain model, not just the same `j`).
    ///
    /// Errors: [`IsogenyEvalError::NoCertificationPoint`] if no point of
    /// order `> 4ℓ²` is found among the first `2^16` x-coordinates (only
    /// for very small `p`); [`IsogenyEvalError::NoDual`] if no candidate
    /// passes the certificate (impossible for a correct `φ`).
    pub fn dual<C: OpCounter>(&self, f: &PrimeField<C>) -> Result<SmallIsogeny, IsogenyEvalError> {
        let q = self
            .certification_point(f)
            .ok_or(IsogenyEvalError::NoCertificationPoint)?;
        let target = self.domain.mul(f, Some(q), self.ell as u128);
        let image = self.eval(f, Some(q));
        let j_domain = self.domain.j_invariant(f);
        for y0 in kernel_x_coordinates(f, &self.codomain, self.ell)? {
            let Ok(back) = SmallIsogeny::new(f, &self.codomain, self.ell, y0) else {
                continue;
            };
            if back.codomain.j_invariant(f) != j_domain {
                continue;
            }
            for u in isomorphisms(f, &back.codomain, &self.domain) {
                let u2 = f.sqr(u);
                let cand = SmallIsogeny {
                    iso: Some((u2, f.mul(u2, u))),
                    codomain: self.domain,
                    ..back
                };
                if cand.eval(f, image) == target {
                    return Ok(cand);
                }
            }
        }
        Err(IsogenyEvalError::NoDual)
    }
}
