//! Certified same-field genus-3 covers of ordinary binary elliptic curves.
//!
//! This module deliberately implements one construction family only.  Over
//! `k = GF(2^m)`, let
//!
//! ```text
//! E : y^2 + x*y = x^3 + a*x^2 + b,        b != 0,
//! d = sqrt(b),
//! H : v^2 + u^2*v = u^7 + a*u^4 + d*u.
//! ```
//!
//! Then
//!
//! ```text
//! phi : H -> E,       (u, v) |-> (u^3, u*v + d)
//! ```
//!
//! is a separable degree-three morphism over `k`.  The smooth projective
//! model of `H` has genus three and one rational point at infinity.  The map
//! is tamely ramified with index three at infinity and at `(u,v)=(0,0)`;
//! the latter maps to the rational two-torsion point `(0,d)` on `E`.
//!
//! The implementation constructs and validates the equations, maps rational
//! points, constructs `phi^*([P]-[O])` in reduced Mumford form, and pushes a
//! rational point divisor `[Q]-[infinity]` forward.  A complete arbitrary-class
//! pushforward is also available for the exact exhaustive `GF(8)` fixture
//! `a=1`, `b=3`, using a checked 486-class presentation; every other model
//! remains explicitly unsupported.  This module contains no discrete-log
//! solver or performance claim.

use super::{
    hyperelliptic::HyperellipticError, BinaryPoint, F2mElement, F2mPoly, HyperellipticCurve,
    IrreduciblePoly, MumfordDivisor,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::sync::OnceLock;

/// Construction and transfer failures.  Unsupported operations are distinct
/// from invalid mathematical input.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CoverTransferError {
    InvalidField(&'static str),
    FieldWidthMismatch {
        object: &'static str,
        expected: u32,
        actual: u32,
    },
    SingularTarget,
    PointNotOnTarget,
    PointNotOnCover,
    InvalidSubgroupOrder,
    DivisorPolynomialTooLong {
        polynomial: &'static str,
        maximum_coefficients: usize,
        actual_coefficients: usize,
    },
    CertificateFailure(&'static str),
    Hyperelliptic(HyperellipticError),
    Unsupported(UnsupportedCoverOperation),
}

impl fmt::Display for CoverTransferError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(reason) => write!(f, "invalid binary field: {reason}"),
            Self::FieldWidthMismatch {
                object,
                expected,
                actual,
            } => write!(f, "{object} has field width {actual}, expected {expected}"),
            Self::SingularTarget => write!(f, "ordinary binary target requires b != 0"),
            Self::PointNotOnTarget => write!(f, "point is not on the target elliptic curve"),
            Self::PointNotOnCover => write!(f, "point is not on the source cover"),
            Self::InvalidSubgroupOrder => write!(f, "subgroup order must be positive"),
            Self::DivisorPolynomialTooLong {
                polynomial,
                maximum_coefficients,
                actual_coefficients,
            } => write!(
                f,
                "fixed-fixture divisor polynomial {polynomial} has {actual_coefficients} coefficients, maximum {maximum_coefficients}"
            ),
            Self::CertificateFailure(reason) => write!(f, "cover certificate failed: {reason}"),
            Self::Hyperelliptic(error) => write!(f, "hyperelliptic validation failed: {error}"),
            Self::Unsupported(op) => write!(f, "unsupported cover operation: {op}"),
        }
    }
}

impl std::error::Error for CoverTransferError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Hyperelliptic(error) => Some(error),
            _ => None,
        }
    }
}

impl From<HyperellipticError> for CoverTransferError {
    fn from(error: HyperellipticError) -> Self {
        Self::Hyperelliptic(error)
    }
}

/// Operations intentionally outside the certified implementation boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum UnsupportedCoverOperation {
    /// The exhaustive arbitrary-class implementation is certified only for
    /// one fixed `GF(8)` fixture.  Other models still require a general
    /// polynomial norm/resultant implementation, including non-split support.
    ArbitraryJacobianPushforward,
}

impl fmt::Display for UnsupportedCoverOperation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ArbitraryJacobianPushforward => {
                write!(f, "pushforward of an arbitrary Jacobian class")
            }
        }
    }
}

/// A rational point on the genus-three source, including its unique point at
/// infinity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BinaryCoverPoint {
    Infinity,
    Affine { u: F2mElement, v: F2mElement },
}

/// Replayable fixed-family geometry facts checked during construction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BinaryCoverCertificate {
    pub construction: &'static str,
    pub field_degree: u32,
    pub modulus_irreducible: bool,
    pub target_nonsingular: bool,
    pub d_square_equals_b: bool,
    pub defining_equation_substitution: bool,
    pub source_model_validated: bool,
    pub source_genus: u32,
    pub map_degree: u32,
    pub map_separable: bool,
    pub unique_rational_infinity: bool,
    pub affine_smooth_at_u_zero: bool,
    pub infinity_smooth: bool,
    pub infinity_maps_to_target_infinity: bool,
    pub infinity_ramification_index: u32,
    pub exceptional_maps_to_two_torsion: bool,
    pub exceptional_ramification_index: u32,
    pub riemann_hurwitz_ramification_total: u32,
    pub artin_schreier_pole_orders: [u32; 2],
}

/// The three possible shapes of a point-class pullback.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PullbackFiberKind {
    InfinityRamified { ramification_index: u32 },
    ExceptionalTwoTorsionRamified { ramification_index: u32 },
    GenericSeparable { geometric_degree: u32 },
}

/// A computable witness for `phi_* phi^* = [3]` on a rational elliptic
/// point class.  This specialized certificate does not claim a general
/// Jacobian pushforward implementation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PullPushCompositionCertificate {
    pub input: BinaryPoint,
    pub pullback: MumfordDivisor,
    pub fiber_kind: PullbackFiberKind,
    pub formal_pushforward: BinaryPoint,
    pub three_times_input: BinaryPoint,
    pub identity_holds: bool,
}

/// Conditional subgroup-preservation consequence of `phi_* phi^* = [3]`.
///
/// This type does not certify that the caller's integer is the order of a
/// particular subgroup or generator.  `InjectiveIfExactSubgroupOrder` applies
/// only after that exact-order premise has been established independently.
/// If `3` divides the order, the composition identity alone does not settle
/// the kernel intersection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SubgroupTransferStatus {
    InjectiveIfExactSubgroupOrder {
        order: BigUint,
        degree_inverse_mod_order: BigUint,
    },
    KernelIntersectionUnresolved {
        order: BigUint,
        common_factor: u32,
    },
}

/// The checked fixed-family cover and its exact field/model data.
#[derive(Clone, Debug)]
pub struct OrdinaryBinaryCover {
    m: u32,
    irreducible: IrreduciblePoly,
    a: F2mElement,
    b: F2mElement,
    d: F2mElement,
    source: HyperellipticCurve,
    certificate: BinaryCoverCertificate,
}

type FixedF8DivisorKey = (Vec<BigUint>, Vec<BigUint>);

#[derive(Clone, Debug)]
struct FixedF8NormTable {
    images: HashMap<FixedF8DivisorKey, BinaryPoint>,
}

static FIXED_F8_NORM_TABLE: OnceLock<Result<FixedF8NormTable, CoverTransferError>> =
    OnceLock::new();

impl OrdinaryBinaryCover {
    /// Construct and validate the fixed ordinary-binary cover family.
    ///
    /// The polynomial basis is accepted only after an exact Rabin
    /// irreducibility check.  Field elements must already use that field
    /// width; unlike `F2mElement::from_biguint`, this API never truncates an
    /// input coefficient.
    pub fn new(
        m: u32,
        irreducible: IrreduciblePoly,
        a: F2mElement,
        b: F2mElement,
    ) -> Result<Self, CoverTransferError> {
        validate_field(m, &irreducible)?;
        check_width("a", &a, m)?;
        check_width("b", &b, m)?;
        if b.is_zero() {
            return Err(CoverTransferError::SingularTarget);
        }

        // The squaring map is an automorphism of GF(2^m), and its inverse is
        // the (m-1)-fold Frobenius: d = b^(2^(m-1)).
        let d = b.square_k_times(m - 1, &irreducible);
        if d.square(&irreducible) != b {
            return Err(CoverTransferError::CertificateFailure(
                "computed square root does not square to b",
            ));
        }

        let zero = F2mElement::zero(m);
        let one = F2mElement::one(m);
        let h = F2mPoly::from_coeffs(vec![zero.clone(), zero.clone(), one.clone()], m);
        let mut f_coeffs = vec![zero; 8];
        f_coeffs[1] = d.clone();
        f_coeffs[4] = a.clone();
        f_coeffs[7] = one;
        let f = F2mPoly::from_coeffs(f_coeffs, m);
        let source = HyperellipticCurve::try_new(m, irreducible.clone(), h, f, 3)?;

        if !verify_substitution(&source, &a, &b, &d) {
            return Err(CoverTransferError::CertificateFailure(
                "defining-equation substitution failed",
            ));
        }

        // For this exact equation, h'=0 and f'=u^6+d.  The only affine
        // zero of h is u=0, where f'=d != 0.  At infinity, t=1/u and
        // w=v/u^4 give w^2+t^2*w=t+a*t^4+d*t^7, whose t-derivative at
        // (0,0) is one.  The two AS poles have order three.
        let certificate = BinaryCoverCertificate {
            construction: "binary_cubic_pullback_v1",
            field_degree: m,
            modulus_irreducible: true,
            target_nonsingular: true,
            d_square_equals_b: true,
            defining_equation_substitution: true,
            source_model_validated: true,
            source_genus: 3,
            map_degree: 3,
            map_separable: true,
            unique_rational_infinity: true,
            affine_smooth_at_u_zero: !d.is_zero(),
            infinity_smooth: true,
            infinity_maps_to_target_infinity: true,
            infinity_ramification_index: 3,
            exceptional_maps_to_two_torsion: true,
            exceptional_ramification_index: 3,
            riemann_hurwitz_ramification_total: 4,
            artin_schreier_pole_orders: [3, 3],
        };
        if !certificate.affine_smooth_at_u_zero {
            return Err(CoverTransferError::CertificateFailure(
                "source is singular at u=0",
            ));
        }

        let cover = Self {
            m,
            irreducible,
            a,
            b,
            d,
            source,
            certificate,
        };
        cover.verify_certificate()?;
        Ok(cover)
    }

    pub fn field_degree(&self) -> u32 {
        self.m
    }

    pub fn irreducible_polynomial(&self) -> &IrreduciblePoly {
        &self.irreducible
    }

    pub fn target_a(&self) -> &F2mElement {
        &self.a
    }

    pub fn target_b(&self) -> &F2mElement {
        &self.b
    }

    pub fn square_root_b(&self) -> &F2mElement {
        &self.d
    }

    pub fn source_curve(&self) -> &HyperellipticCurve {
        &self.source
    }

    pub fn certificate(&self) -> &BinaryCoverCertificate {
        &self.certificate
    }

    /// Recompute every executable part of the fixed-family certificate.
    ///
    /// The infinity and ramification checks are deductions from the exact
    /// equation shape checked here: the infinity chart is
    /// `w^2+t^2*w=t+a*t^4+d*t^7`, while the Artin--Schreier form has pole
    /// order three at `u=0` and `u=infinity`.
    pub fn verify_certificate(&self) -> Result<(), CoverTransferError> {
        validate_field(self.m, &self.irreducible)?;
        check_width("a", &self.a, self.m)?;
        check_width("b", &self.b, self.m)?;
        check_width("d", &self.d, self.m)?;
        if self.b.is_zero() {
            return Err(CoverTransferError::SingularTarget);
        }
        if self.d.square(&self.irreducible) != self.b {
            return Err(CoverTransferError::CertificateFailure(
                "d^2 does not equal b",
            ));
        }
        self.source.validate()?;
        if self.source.m != self.m
            || self.source.irr.degree != self.irreducible.degree
            || self.source.irr.low_terms != self.irreducible.low_terms
            || self.source.genus != 3
        {
            return Err(CoverTransferError::CertificateFailure(
                "source curve field or genus changed",
            ));
        }

        let zero = F2mElement::zero(self.m);
        let one = F2mElement::one(self.m);
        let expected_h =
            F2mPoly::from_coeffs(vec![zero.clone(), zero.clone(), one.clone()], self.m);
        let mut expected_f_coefficients = vec![zero; 8];
        expected_f_coefficients[1] = self.d.clone();
        expected_f_coefficients[4] = self.a.clone();
        expected_f_coefficients[7] = one;
        let expected_f = F2mPoly::from_coeffs(expected_f_coefficients, self.m);
        if self.source.h != expected_h || self.source.f != expected_f {
            return Err(CoverTransferError::CertificateFailure(
                "source equation is not the certified fixed family",
            ));
        }
        if !verify_substitution(&self.source, &self.a, &self.b, &self.d) {
            return Err(CoverTransferError::CertificateFailure(
                "defining-equation substitution failed",
            ));
        }

        let exceptional_source = BinaryCoverPoint::Affine {
            u: F2mElement::zero(self.m),
            v: F2mElement::zero(self.m),
        };
        let exceptional_target = BinaryPoint::Affine {
            x: F2mElement::zero(self.m),
            y: self.d.clone(),
        };
        if self.map_point(&exceptional_source)? != exceptional_target
            || self.target_add(&exceptional_target, &exceptional_target)? != BinaryPoint::Infinity
            || self.map_point(&BinaryCoverPoint::Infinity)? != BinaryPoint::Infinity
        {
            return Err(CoverTransferError::CertificateFailure(
                "branch-point map certificate failed",
            ));
        }

        let cert = &self.certificate;
        let riemann_hurwitz_left = 2 * cert.source_genus - 2;
        let riemann_hurwitz_right =
            (cert.infinity_ramification_index - 1) + (cert.exceptional_ramification_index - 1);
        if !cert.modulus_irreducible
            || !cert.target_nonsingular
            || !cert.d_square_equals_b
            || !cert.defining_equation_substitution
            || !cert.source_model_validated
            || cert.field_degree != self.m
            || cert.source_genus != 3
            || cert.map_degree != 3
            || !cert.map_separable
            || !cert.unique_rational_infinity
            || !cert.affine_smooth_at_u_zero
            || !cert.infinity_smooth
            || !cert.infinity_maps_to_target_infinity
            || cert.infinity_ramification_index != 3
            || !cert.exceptional_maps_to_two_torsion
            || cert.exceptional_ramification_index != 3
            || cert.artin_schreier_pole_orders != [3, 3]
            || cert.riemann_hurwitz_ramification_total != riemann_hurwitz_right
            || riemann_hurwitz_left != riemann_hurwitz_right
        {
            return Err(CoverTransferError::CertificateFailure(
                "fixed-family geometry certificate is inconsistent",
            ));
        }
        Ok(())
    }

    /// Check the normalized ordinary target equation.
    pub fn target_is_on_curve(&self, point: &BinaryPoint) -> bool {
        match point {
            BinaryPoint::Infinity => true,
            BinaryPoint::Affine { x, y } => {
                if x.m_value() != self.m || y.m_value() != self.m {
                    return false;
                }
                let lhs = y
                    .square(&self.irreducible)
                    .add(&x.mul(y, &self.irreducible));
                let x2 = x.square(&self.irreducible);
                let rhs = x2
                    .mul(x, &self.irreducible)
                    .add(&self.a.mul(&x2, &self.irreducible))
                    .add(&self.b);
                lhs == rhs
            }
        }
    }

    /// Check the affine source equation, with infinity handled explicitly.
    pub fn source_is_on_curve(&self, point: &BinaryCoverPoint) -> bool {
        match point {
            BinaryCoverPoint::Infinity => true,
            BinaryCoverPoint::Affine { u, v } => {
                if u.m_value() != self.m || v.m_value() != self.m {
                    return false;
                }
                self.source.is_on_curve(u, v)
            }
        }
    }

    /// Evaluate `phi` on a rational source point.
    pub fn map_point(&self, point: &BinaryCoverPoint) -> Result<BinaryPoint, CoverTransferError> {
        if !self.source_is_on_curve(point) {
            return Err(CoverTransferError::PointNotOnCover);
        }
        let image = match point {
            BinaryCoverPoint::Infinity => BinaryPoint::Infinity,
            BinaryCoverPoint::Affine { u, v } => {
                let u2 = u.square(&self.irreducible);
                let x = u2.mul(u, &self.irreducible);
                let y = u.mul(v, &self.irreducible).add(&self.d);
                BinaryPoint::Affine { x, y }
            }
        };
        if !self.target_is_on_curve(&image) {
            return Err(CoverTransferError::CertificateFailure(
                "point map failed target equation",
            ));
        }
        Ok(image)
    }

    /// Push `[Q]-[infinity_H]` forward for one rational source point.
    /// The result is represented by the corresponding point of `E`, with
    /// `Infinity` denoting the identity class.
    pub fn pushforward_rational_point_divisor(
        &self,
        point: &BinaryCoverPoint,
    ) -> Result<BinaryPoint, CoverTransferError> {
        self.map_point(point)
    }

    /// Pull `[P]-[O_E]` back and reduce it in `Jac(H)`.
    ///
    /// For `x != 0`, the degree-three fiber has
    /// `U(t)=t^3+x` and `V(t)=((y+d)/x)t^2`.  At the exceptional
    /// two-torsion point `T=(0,d)`, the scheme fiber is `3(0,0)` and its
    /// reduced class is `(U,V)=(t,0)`.
    pub fn pullback_point_class(
        &self,
        point: &BinaryPoint,
    ) -> Result<MumfordDivisor, CoverTransferError> {
        if !self.target_is_on_curve(point) {
            return Err(CoverTransferError::PointNotOnTarget);
        }
        let divisor = match point {
            BinaryPoint::Infinity => MumfordDivisor::identity(self.m),
            BinaryPoint::Affine { x, y } if x.is_zero() => {
                if *y != self.d {
                    return Err(CoverTransferError::CertificateFailure(
                        "the unique point above x=0 is not (0,d)",
                    ));
                }
                MumfordDivisor::try_new(&self.source, F2mPoly::x(self.m), F2mPoly::zero(self.m))?
            }
            BinaryPoint::Affine { x, y } => {
                let x_inv = x.flt_inverse(&self.irreducible).ok_or(
                    CoverTransferError::CertificateFailure("nonzero x is not invertible"),
                )?;
                let c = y.add(&self.d).mul(&x_inv, &self.irreducible);
                let zero = F2mElement::zero(self.m);
                let one = F2mElement::one(self.m);
                MumfordDivisor::try_new(
                    &self.source,
                    F2mPoly::from_coeffs(vec![x.clone(), zero.clone(), zero.clone(), one], self.m),
                    F2mPoly::from_coeffs(vec![zero.clone(), zero, c], self.m),
                )?
            }
        };
        if !divisor.is_reduced(&self.source) {
            return Err(CoverTransferError::CertificateFailure(
                "pullback did not produce a reduced valid Mumford divisor",
            ));
        }
        Ok(divisor)
    }

    /// Return the fiber type used by the point-class pullback certificate.
    pub fn pullback_fiber_kind(
        &self,
        point: &BinaryPoint,
    ) -> Result<PullbackFiberKind, CoverTransferError> {
        if !self.target_is_on_curve(point) {
            return Err(CoverTransferError::PointNotOnTarget);
        }
        Ok(match point {
            BinaryPoint::Infinity => PullbackFiberKind::InfinityRamified {
                ramification_index: 3,
            },
            BinaryPoint::Affine { x, .. } if x.is_zero() => {
                PullbackFiberKind::ExceptionalTwoTorsionRamified {
                    ramification_index: 3,
                }
            }
            BinaryPoint::Affine { .. } => PullbackFiberKind::GenericSeparable {
                geometric_degree: 3,
            },
        })
    }

    /// Build a specialized, replayable certificate for
    /// `phi_* phi^*([P]-[O]) = [3]([P]-[O])`.
    pub fn pull_push_composition_certificate(
        &self,
        point: &BinaryPoint,
    ) -> Result<PullPushCompositionCertificate, CoverTransferError> {
        let pullback = self.pullback_point_class(point)?;
        let fiber_kind = self.pullback_fiber_kind(point)?;
        let formal_pushforward =
            self.formal_pushforward_of_point_pullback(point, &pullback, &fiber_kind)?;
        let three_times_input = self.target_scalar_mul_small(point, 3)?;
        let identity_holds = formal_pushforward == three_times_input;
        if !identity_holds {
            return Err(CoverTransferError::CertificateFailure(
                "push-pull composition did not equal multiplication by three",
            ));
        }
        Ok(PullPushCompositionCertificate {
            input: point.clone(),
            pullback,
            fiber_kind,
            formal_pushforward,
            three_times_input,
            identity_holds,
        })
    }

    /// Push the typed pullback witness forward without invoking a general
    /// Jacobian norm.  This is deliberately private: it is valid only for
    /// the three point fibers constructed by `pullback_point_class`.
    fn formal_pushforward_of_point_pullback(
        &self,
        point: &BinaryPoint,
        pullback: &MumfordDivisor,
        fiber_kind: &PullbackFiberKind,
    ) -> Result<BinaryPoint, CoverTransferError> {
        match (point, fiber_kind) {
            (
                BinaryPoint::Infinity,
                PullbackFiberKind::InfinityRamified {
                    ramification_index: 3,
                },
            ) => {
                if !pullback.is_identity() {
                    return Err(CoverTransferError::CertificateFailure(
                        "infinity pullback is not the identity class",
                    ));
                }
                self.pushforward_rational_point_divisor(&BinaryCoverPoint::Infinity)
            }
            (
                BinaryPoint::Affine { x, y },
                PullbackFiberKind::ExceptionalTwoTorsionRamified {
                    ramification_index: 3,
                },
            ) if x.is_zero() && *y == self.d => {
                let expected = MumfordDivisor::try_new(
                    &self.source,
                    F2mPoly::x(self.m),
                    F2mPoly::zero(self.m),
                )?;
                if pullback != &expected {
                    return Err(CoverTransferError::CertificateFailure(
                        "exceptional reduced pullback witness changed",
                    ));
                }
                if !expected.try_double(&self.source)?.is_identity() {
                    return Err(CoverTransferError::CertificateFailure(
                        "exceptional source-point class is not two-torsion",
                    ));
                }
                let image = self.pushforward_rational_point_divisor(&BinaryCoverPoint::Affine {
                    u: F2mElement::zero(self.m),
                    v: F2mElement::zero(self.m),
                })?;
                // Scheme-theoretically, phi^*(T-O)=3D for
                // D=[(0,0)]-[infinity].  The checked equality 2D=0 makes
                // the returned reduced pair represent 3D=D.  We push this
                // one rational-point class only; this is not a general
                // Jacobian norm.  Its image T=(0,d) also has order two, so
                // phi_*(3D)=phi_*(D)=T=[3]T.
                if self.target_scalar_mul_small(&image, 3)? != image {
                    return Err(CoverTransferError::CertificateFailure(
                        "exceptional image is not fixed by multiplication by three",
                    ));
                }
                Ok(image)
            }
            (
                BinaryPoint::Affine { x, y },
                PullbackFiberKind::GenericSeparable { geometric_degree },
            ) if !x.is_zero() && *geometric_degree == 3 => {
                let x_inverse = x.flt_inverse(&self.irreducible).ok_or(
                    CoverTransferError::CertificateFailure("generic fiber x is not invertible"),
                )?;
                let coefficient = y.add(&self.d).mul(&x_inverse, &self.irreducible);
                let zero = F2mElement::zero(self.m);
                let expected = MumfordDivisor::try_new(
                    &self.source,
                    F2mPoly::from_coeffs(
                        vec![
                            x.clone(),
                            zero.clone(),
                            zero.clone(),
                            F2mElement::one(self.m),
                        ],
                        self.m,
                    ),
                    F2mPoly::from_coeffs(vec![zero.clone(), zero, coefficient], self.m),
                )?;
                if pullback != &expected {
                    return Err(CoverTransferError::CertificateFailure(
                        "generic pullback witness changed",
                    ));
                }

                // U=t^3+x is separable for x!=0: U'=t^2 has no common
                // root with U.  Its three geometric support points all map
                // to P, so sum their point classes directly.
                let mut image_sum = BinaryPoint::Infinity;
                for _ in 0..*geometric_degree {
                    image_sum = self.target_add(&image_sum, point)?;
                }
                Ok(image_sum)
            }
            _ => Err(CoverTransferError::CertificateFailure(
                "point and pullback fiber kind are inconsistent",
            )),
        }
    }

    /// Derive the subgroup-preservation condition from `phi_* phi^* = [3]`.
    ///
    /// The integer is not bound to a subgroup or generator by this method.
    /// Therefore the injective result is conditional on the caller separately
    /// proving that `order` is the exact order of the subgroup in question.
    pub fn subgroup_transfer_status(
        &self,
        order: &BigUint,
    ) -> Result<SubgroupTransferStatus, CoverTransferError> {
        if order.is_zero() {
            return Err(CoverTransferError::InvalidSubgroupOrder);
        }
        let three = BigUint::from(3u32);
        let residue = order % &three;
        if residue.is_zero() {
            return Ok(SubgroupTransferStatus::KernelIntersectionUnresolved {
                order: order.clone(),
                common_factor: 3,
            });
        }
        let inverse = if residue == BigUint::one() {
            (order * 2u32 + 1u32) / 3u32
        } else {
            (order + 1u32) / 3u32
        };
        debug_assert_eq!((&inverse * 3u32) % order, BigUint::one() % order);
        Ok(SubgroupTransferStatus::InjectiveIfExactSubgroupOrder {
            order: order.clone(),
            degree_inverse_mod_order: inverse,
        })
    }

    /// Push an arbitrary reduced Jacobian class to the target elliptic curve.
    ///
    /// This operation is deliberately bounded to the exact exhaustive fixture
    ///
    /// ```text
    /// GF(8) = GF(2)[z]/(z^3+z+1),  a=1, b=3,
    /// H: v^2+u^2*v=u^7+u^4+7u.
    /// ```
    ///
    /// Its 486 rational Jacobian classes have the unique presentation
    /// `i*A+j*B`, `0<=i<9`, `0<=j<54`, where
    /// `A=[(2,3)-infinity]` and `B=[(3,1)-infinity]`.  The cached table is
    /// constructed with the checked Jacobian law and assigns each class the
    /// corresponding sum of the two directly mapped rational point classes.
    /// This handles non-split support in this finite fixture without claiming
    /// a general resultant implementation.  Every other cover returns a
    /// structured unsupported result.
    pub fn pushforward_arbitrary_jacobian_class(
        &self,
        divisor: &MumfordDivisor,
    ) -> Result<BinaryPoint, CoverTransferError> {
        if !self.is_fixed_f8_norm_fixture() {
            return Err(CoverTransferError::Unsupported(
                UnsupportedCoverOperation::ArbitraryJacobianPushforward,
            ));
        }
        // A reduced genus-three Mumford representative has deg(u) <= 3 and
        // deg(v) < deg(u), hence at most four and three coefficients.  Check
        // the public vector lengths before the full curve/divisor validation,
        // whose canonicality and field-width checks otherwise scan every
        // supplied coefficient.
        if divisor.u.coeffs.len() > 4 {
            return Err(CoverTransferError::DivisorPolynomialTooLong {
                polynomial: "u",
                maximum_coefficients: 4,
                actual_coefficients: divisor.u.coeffs.len(),
            });
        }
        if divisor.v.coeffs.len() > 3 {
            return Err(CoverTransferError::DivisorPolynomialTooLong {
                polynomial: "v",
                maximum_coefficients: 3,
                actual_coefficients: divisor.v.coeffs.len(),
            });
        }
        divisor.validate(&self.source)?;
        let table = FIXED_F8_NORM_TABLE
            .get_or_init(|| self.build_fixed_f8_norm_table())
            .as_ref()
            .map_err(Clone::clone)?;
        table
            .images
            .get(&fixed_f8_divisor_key(divisor))
            .cloned()
            .ok_or(CoverTransferError::CertificateFailure(
                "validated fixed-fixture divisor is absent from the 486-class norm table",
            ))
    }

    fn is_fixed_f8_norm_fixture(&self) -> bool {
        self.m == 3
            && self.irreducible.degree == 3
            && self.irreducible.low_terms.as_slice() == [0u32, 1u32]
            && self.a == F2mElement::one(3)
            && self.b == F2mElement::from_biguint(&BigUint::from(3u32), 3)
            && self.d == F2mElement::from_biguint(&BigUint::from(7u32), 3)
    }

    fn build_fixed_f8_norm_table(&self) -> Result<FixedF8NormTable, CoverTransferError> {
        if !self.is_fixed_f8_norm_fixture() {
            return Err(CoverTransferError::Unsupported(
                UnsupportedCoverOperation::ArbitraryJacobianPushforward,
            ));
        }
        self.verify_certificate()?;

        let a_u = F2mElement::from_biguint(&BigUint::from(2u32), 3);
        let a_v = F2mElement::from_biguint(&BigUint::from(3u32), 3);
        let b_u = F2mElement::from_biguint(&BigUint::from(3u32), 3);
        let b_v = F2mElement::one(3);
        let a = MumfordDivisor::try_from_point(&self.source, &a_u, &a_v)?;
        let b = MumfordDivisor::try_from_point(&self.source, &b_u, &b_v)?;
        let a_point = BinaryCoverPoint::Affine { u: a_u, v: a_v };
        let b_point = BinaryCoverPoint::Affine { u: b_u, v: b_v };
        let image_a = self.pushforward_rational_point_divisor(&a_point)?;
        let image_b = self.pushforward_rational_point_divisor(&b_point)?;
        let expected_image_a = BinaryPoint::Affine {
            x: F2mElement::from_biguint(&BigUint::from(3u32), 3),
            y: F2mElement::one(3),
        };
        let expected_image_b = BinaryPoint::Affine {
            x: F2mElement::from_biguint(&BigUint::from(4u32), 3),
            y: F2mElement::from_biguint(&BigUint::from(4u32), 3),
        };

        if !a
            .try_scalar_mul(&BigUint::from(27u32), &self.source)?
            .is_identity()
            || a.try_scalar_mul(&BigUint::from(9u32), &self.source)?
                .is_identity()
        {
            return Err(CoverTransferError::CertificateFailure(
                "fixed F8 norm generator A does not have order 27",
            ));
        }
        if !b
            .try_scalar_mul(&BigUint::from(54u32), &self.source)?
            .is_identity()
            || b.try_scalar_mul(&BigUint::from(27u32), &self.source)?
                .is_identity()
            || b.try_scalar_mul(&BigUint::from(18u32), &self.source)?
                .is_identity()
        {
            return Err(CoverTransferError::CertificateFailure(
                "fixed F8 norm generator B does not have order 54",
            ));
        }
        if a.try_scalar_mul(&BigUint::from(9u32), &self.source)?
            != b.try_scalar_mul(&BigUint::from(18u32), &self.source)?
        {
            return Err(CoverTransferError::CertificateFailure(
                "fixed F8 norm generators do not satisfy 9A=18B",
            ));
        }
        if image_a != expected_image_a
            || image_b != expected_image_b
            || image_a != self.target_scalar_mul_small(&image_b, 2)?
            || self.target_scalar_mul_small(&image_b, 6)? != BinaryPoint::Infinity
            || self.target_scalar_mul_small(&image_b, 3)? == BinaryPoint::Infinity
            || self.target_scalar_mul_small(&image_b, 2)? == BinaryPoint::Infinity
        {
            return Err(CoverTransferError::CertificateFailure(
                "fixed F8 point images do not satisfy phi(A)=2Q with Q of order six",
            ));
        }

        let a_multiples = self.fixed_f8_jacobian_multiples(&a, 9)?;
        let b_multiples = self.fixed_f8_jacobian_multiples(&b, 54)?;
        let image_a_multiples = self.fixed_f8_target_multiples(&image_a, 9)?;
        let image_b_multiples = self.fixed_f8_target_multiples(&image_b, 54)?;
        let mut images = HashMap::with_capacity(486);
        for (a_multiple, image_a_multiple) in a_multiples.iter().zip(image_a_multiples.iter()) {
            for (b_multiple, image_b_multiple) in b_multiples.iter().zip(image_b_multiples.iter()) {
                let divisor = a_multiple.try_add(b_multiple, &self.source)?;
                let image = self.target_add(image_a_multiple, image_b_multiple)?;
                let key = fixed_f8_divisor_key(&divisor);
                if let Some(previous) = images.insert(key, image.clone()) {
                    if previous != image {
                        return Err(CoverTransferError::CertificateFailure(
                            "fixed F8 norm presentation assigns inconsistent images",
                        ));
                    }
                }
            }
        }
        if images.len() != 486 {
            return Err(CoverTransferError::CertificateFailure(
                "fixed F8 norm presentation does not contain 486 unique classes",
            ));
        }
        Ok(FixedF8NormTable { images })
    }

    fn fixed_f8_jacobian_multiples(
        &self,
        base: &MumfordDivisor,
        count: usize,
    ) -> Result<Vec<MumfordDivisor>, CoverTransferError> {
        let mut multiples = Vec::with_capacity(count);
        let mut current = MumfordDivisor::identity(self.m);
        for _ in 0..count {
            multiples.push(current.clone());
            current = current.try_add(base, &self.source)?;
        }
        Ok(multiples)
    }

    fn fixed_f8_target_multiples(
        &self,
        base: &BinaryPoint,
        count: usize,
    ) -> Result<Vec<BinaryPoint>, CoverTransferError> {
        let mut multiples = Vec::with_capacity(count);
        let mut current = BinaryPoint::Infinity;
        for _ in 0..count {
            multiples.push(current.clone());
            current = self.target_add(&current, base)?;
        }
        Ok(multiples)
    }

    fn target_scalar_mul_small(
        &self,
        point: &BinaryPoint,
        scalar: u32,
    ) -> Result<BinaryPoint, CoverTransferError> {
        if !self.target_is_on_curve(point) {
            return Err(CoverTransferError::PointNotOnTarget);
        }
        let mut acc = BinaryPoint::Infinity;
        let mut base = point.clone();
        let mut k = scalar;
        while k != 0 {
            if k & 1 == 1 {
                acc = self.target_add(&acc, &base)?;
            }
            k >>= 1;
            if k != 0 {
                base = self.target_add(&base, &base)?;
            }
        }
        Ok(acc)
    }

    fn target_add(
        &self,
        left: &BinaryPoint,
        right: &BinaryPoint,
    ) -> Result<BinaryPoint, CoverTransferError> {
        if !self.target_is_on_curve(left) || !self.target_is_on_curve(right) {
            return Err(CoverTransferError::PointNotOnTarget);
        }
        let result = match (left, right) {
            (BinaryPoint::Infinity, p) | (p, BinaryPoint::Infinity) => p.clone(),
            (BinaryPoint::Affine { x: x1, y: y1 }, BinaryPoint::Affine { x: x2, y: y2 }) => {
                if x1 == x2 {
                    if y1.add(y2) == *x1 || x1.is_zero() {
                        BinaryPoint::Infinity
                    } else {
                        let inv = x1.flt_inverse(&self.irreducible).ok_or(
                            CoverTransferError::CertificateFailure(
                                "doubling denominator is not invertible",
                            ),
                        )?;
                        let lambda = x1.add(&y1.mul(&inv, &self.irreducible));
                        let x3 = lambda.square(&self.irreducible).add(&lambda).add(&self.a);
                        let y3 = x1.square(&self.irreducible).add(
                            &lambda
                                .add(&F2mElement::one(self.m))
                                .mul(&x3, &self.irreducible),
                        );
                        BinaryPoint::Affine { x: x3, y: y3 }
                    }
                } else {
                    let denominator = x1.add(x2);
                    let inv = denominator.flt_inverse(&self.irreducible).ok_or(
                        CoverTransferError::CertificateFailure(
                            "addition denominator is not invertible",
                        ),
                    )?;
                    let lambda = y1.add(y2).mul(&inv, &self.irreducible);
                    let x3 = lambda
                        .square(&self.irreducible)
                        .add(&lambda)
                        .add(x1)
                        .add(x2)
                        .add(&self.a);
                    let y3 = lambda.mul(&x1.add(&x3), &self.irreducible).add(&x3).add(y1);
                    BinaryPoint::Affine { x: x3, y: y3 }
                }
            }
        };
        if !self.target_is_on_curve(&result) {
            return Err(CoverTransferError::CertificateFailure(
                "target group operation left the curve",
            ));
        }
        Ok(result)
    }
}

fn fixed_f8_divisor_key(divisor: &MumfordDivisor) -> FixedF8DivisorKey {
    (
        divisor
            .u
            .coeffs
            .iter()
            .map(F2mElement::to_biguint)
            .collect(),
        divisor
            .v
            .coeffs
            .iter()
            .map(F2mElement::to_biguint)
            .collect(),
    )
}

fn check_width(
    object: &'static str,
    value: &F2mElement,
    expected: u32,
) -> Result<(), CoverTransferError> {
    let actual = value.m_value();
    if actual != expected {
        return Err(CoverTransferError::FieldWidthMismatch {
            object,
            expected,
            actual,
        });
    }
    Ok(())
}

fn validate_field(m: u32, polynomial: &IrreduciblePoly) -> Result<(), CoverTransferError> {
    if !(1..=4096).contains(&m) {
        return Err(CoverTransferError::InvalidField(
            "degree must lie in 1..=4096",
        ));
    }
    if polynomial.degree != m {
        return Err(CoverTransferError::InvalidField(
            "modulus degree does not equal field width",
        ));
    }
    if polynomial.low_terms.first() != Some(&0) {
        return Err(CoverTransferError::InvalidField(
            "modulus must have nonzero constant term and canonical term order",
        ));
    }
    if polynomial.low_terms.iter().any(|&term| term >= m)
        || polynomial.low_terms.windows(2).any(|w| w[0] >= w[1])
    {
        return Err(CoverTransferError::InvalidField(
            "modulus low terms must be unique, increasing, and below m",
        ));
    }
    if !rabin_irreducible(polynomial) {
        return Err(CoverTransferError::InvalidField(
            "modulus is reducible over GF(2)",
        ));
    }
    Ok(())
}

fn modulus_bits(polynomial: &IrreduciblePoly) -> BigUint {
    let mut modulus = BigUint::one() << polynomial.degree as usize;
    for &term in &polynomial.low_terms {
        modulus |= BigUint::one() << term as usize;
    }
    modulus
}

fn bit_remainder(mut value: BigUint, modulus: &BigUint) -> BigUint {
    while !value.is_zero() && value.bits() >= modulus.bits() {
        let shift = (value.bits() - modulus.bits()) as usize;
        value ^= modulus << shift;
    }
    value
}

fn bit_gcd(mut left: BigUint, mut right: BigUint) -> BigUint {
    while !right.is_zero() {
        let remainder = bit_remainder(left, &right);
        left = right;
        right = remainder;
    }
    left
}

/// Rabin's exact irreducibility criterion over GF(2).
fn rabin_irreducible(polynomial: &IrreduciblePoly) -> bool {
    let m = polynomial.degree;
    let modulus = modulus_bits(polynomial);
    let mut remaining = m;
    let mut checkpoints = HashSet::new();
    let mut prime = 2u32;
    while prime * prime <= remaining {
        if remaining.is_multiple_of(prime) {
            checkpoints.insert(m / prime);
            while remaining.is_multiple_of(prime) {
                remaining /= prime;
            }
        }
        prime += 1;
    }
    if remaining > 1 {
        checkpoints.insert(m / remaining);
    }

    let z = bit_remainder(BigUint::from(2u32), &modulus);
    let mut frobenius = F2mElement::from_biguint(&z, m);
    for i in 1..=m {
        frobenius = frobenius.square(polynomial);
        if checkpoints.contains(&i)
            && bit_gcd(frobenius.to_biguint() ^ &z, modulus.clone()) != BigUint::one()
        {
            return false;
        }
    }
    frobenius.to_biguint() == z
}

fn verify_substitution(
    source: &HyperellipticCurve,
    a: &F2mElement,
    b: &F2mElement,
    d: &F2mElement,
) -> bool {
    let m = source.m;
    let irr = &source.irr;
    let t = F2mPoly::x(m);
    let t2 = t.mul(&t, irr);
    let x = t2.mul(&t, irr);
    let x2 = x.mul(&x, irr);
    let x3 = x2.mul(&x, irr);
    let aa = t2;

    // y=t*v+d.  Reducing y^2+x*y modulo v^2=h*v+f leaves a coefficient
    // of v and a constant coefficient, both of which must match the target
    // right-hand side.
    let linear = aa.mul(&source.h, irr).add(&x.mul(&t, irr));
    let constant = aa
        .mul(&source.f, irr)
        .add(&F2mPoly::constant(d.square(irr)))
        .add(&x.scalar_mul(d, irr))
        .add(&x3)
        .add(&x2.scalar_mul(a, irr))
        .add(&F2mPoly::constant(b.clone()));
    linear.is_zero() && constant.is_zero()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn f8_modulus() -> IrreduciblePoly {
        // z^3 + z + 1.
        IrreduciblePoly {
            degree: 3,
            low_terms: vec![0, 1],
        }
    }

    fn elt(value: u32) -> F2mElement {
        F2mElement::from_biguint(&BigUint::from(value), 3)
    }

    fn cover(a: u32, b: u32) -> OrdinaryBinaryCover {
        OrdinaryBinaryCover::new(3, f8_modulus(), elt(a), elt(b)).unwrap()
    }

    fn coefficients_from_code(mut code: usize, length: usize) -> Vec<F2mElement> {
        (0..length)
            .map(|_| {
                let coefficient = elt((code % 8) as u32);
                code /= 8;
                coefficient
            })
            .collect()
    }

    fn all_fixed_f8_reduced_divisors(curve: &HyperellipticCurve) -> Vec<MumfordDivisor> {
        let mut divisors = vec![MumfordDivisor::identity(3)];
        for u_degree in 1..=3usize {
            let coefficient_vectors = 8usize.pow(u_degree as u32);
            for u_code in 0..coefficient_vectors {
                let mut u_coefficients = coefficients_from_code(u_code, u_degree);
                u_coefficients.push(F2mElement::one(3));
                let u = F2mPoly::from_coeffs(u_coefficients, 3);
                for v_code in 0..coefficient_vectors {
                    let v = F2mPoly::from_coeffs(coefficients_from_code(v_code, u_degree), 3);
                    // These loops already enforce monic U and
                    // deg(V)<deg(U)<=3.  Check the remaining Mumford
                    // invariant directly so exhaustive enumeration does not
                    // revalidate the fixed field for all 266,304 candidates.
                    let equation = v
                        .mul(&v, &curve.irr)
                        .add(&curve.h.mul(&v, &curve.irr))
                        .add(&curve.f);
                    let (_quotient, remainder) = equation.divrem(&u, &curve.irr);
                    if remainder.is_zero() {
                        divisors.push(MumfordDivisor { u: u.clone(), v });
                    }
                }
            }
        }
        divisors
    }

    fn target_point(x: u32, y: u32) -> BinaryPoint {
        BinaryPoint::Affine {
            x: elt(x),
            y: elt(y),
        }
    }

    fn checked_divisor(curve: &HyperellipticCurve, u: &[u32], v: &[u32]) -> MumfordDivisor {
        MumfordDivisor::try_new(
            curve,
            F2mPoly::from_coeffs(u.iter().copied().map(elt).collect(), 3),
            F2mPoly::from_coeffs(v.iter().copied().map(elt).collect(), 3),
        )
        .unwrap()
    }

    #[test]
    fn construction_checks_geometry_and_rejects_invalid_inputs() {
        let c = cover(1, 1);
        let cert = c.certificate();
        assert_eq!(cert.construction, "binary_cubic_pullback_v1");
        assert!(cert.modulus_irreducible);
        assert!(cert.target_nonsingular);
        assert!(cert.d_square_equals_b);
        assert!(cert.defining_equation_substitution);
        assert!(cert.source_model_validated);
        assert_eq!(cert.source_genus, 3);
        assert_eq!(cert.map_degree, 3);
        assert!(cert.map_separable);
        assert!(cert.unique_rational_infinity);
        assert!(cert.affine_smooth_at_u_zero);
        assert!(cert.infinity_smooth);
        assert_eq!(cert.artin_schreier_pole_orders, [3, 3]);
        assert_eq!(cert.riemann_hurwitz_ramification_total, 4);
        assert_eq!(c.field_degree(), 3);
        assert_eq!(c.irreducible_polynomial().degree, 3);
        assert_eq!(c.target_a(), &elt(1));
        assert_eq!(c.target_b(), &elt(1));
        assert_eq!(c.square_root_b(), &elt(1));
        c.verify_certificate().unwrap();

        assert_eq!(
            OrdinaryBinaryCover::new(3, f8_modulus(), elt(0), elt(0)).unwrap_err(),
            CoverTransferError::SingularTarget
        );
        let reducible = IrreduciblePoly {
            degree: 3,
            low_terms: vec![0, 1, 2], // z^3+z^2+z+1=(z+1)^3
        };
        assert!(matches!(
            OrdinaryBinaryCover::new(3, reducible, elt(0), elt(1)),
            Err(CoverTransferError::InvalidField(_))
        ));
        assert!(matches!(
            OrdinaryBinaryCover::new(3, f8_modulus(), F2mElement::one(4), elt(1)),
            Err(CoverTransferError::FieldWidthMismatch { object: "a", .. })
        ));

        let off_target = BinaryPoint::Affine {
            x: elt(0),
            y: elt(0),
        };
        assert_eq!(
            c.pullback_point_class(&off_target).unwrap_err(),
            CoverTransferError::PointNotOnTarget
        );
        let off_source = BinaryCoverPoint::Affine {
            u: elt(0),
            v: elt(1),
        };
        assert_eq!(
            c.map_point(&off_source).unwrap_err(),
            CoverTransferError::PointNotOnCover
        );
    }

    #[test]
    fn infinity_and_exceptional_two_torsion_are_explicit() {
        let c = cover(1, 1);
        assert_eq!(
            c.map_point(&BinaryCoverPoint::Infinity).unwrap(),
            BinaryPoint::Infinity
        );
        let infinity = c
            .pull_push_composition_certificate(&BinaryPoint::Infinity)
            .unwrap();
        assert!(infinity.pullback.is_identity());
        assert_eq!(
            infinity.fiber_kind,
            PullbackFiberKind::InfinityRamified {
                ramification_index: 3
            }
        );
        assert!(infinity.identity_holds);

        let p0 = BinaryCoverPoint::Affine {
            u: elt(0),
            v: elt(0),
        };
        let t = BinaryPoint::Affine {
            x: elt(0),
            y: c.d.clone(),
        };
        assert!(c.source_is_on_curve(&p0));
        assert_eq!(c.map_point(&p0).unwrap(), t);
        let pullback = c.pullback_point_class(&t).unwrap();
        assert_eq!(pullback.u, F2mPoly::x(3));
        assert!(pullback.v.is_zero());
        assert!(pullback.is_reduced(c.source_curve()));
        assert!(pullback.try_double(c.source_curve()).unwrap().is_identity());
        let cert = c.pull_push_composition_certificate(&t).unwrap();
        assert_eq!(
            cert.fiber_kind,
            PullbackFiberKind::ExceptionalTwoTorsionRamified {
                ramification_index: 3
            }
        );
        assert_eq!(cert.formal_pushforward, t);
        assert!(cert.identity_holds);
    }

    #[test]
    fn exhaustive_f8_models_points_maps_and_compositions() {
        let mut target_points = 0usize;
        let mut source_points = 0usize;
        let mut generic_pullbacks = 0usize;
        for a in 0..8 {
            for b in 1..8 {
                let c = cover(a, b);
                for x in 0..8 {
                    for y in 0..8 {
                        let point = BinaryPoint::Affine {
                            x: elt(x),
                            y: elt(y),
                        };
                        if !c.target_is_on_curve(&point) {
                            continue;
                        }
                        target_points += 1;
                        let cert = c.pull_push_composition_certificate(&point).unwrap();
                        assert!(cert.pullback.is_reduced(c.source_curve()));
                        assert!(cert.identity_holds);
                        if x != 0 {
                            assert_eq!(cert.pullback.u.degree(), Some(3));
                            assert!(matches!(
                                cert.fiber_kind,
                                PullbackFiberKind::GenericSeparable {
                                    geometric_degree: 3
                                }
                            ));
                            generic_pullbacks += 1;
                        }
                    }
                }
                let infinity = c
                    .pull_push_composition_certificate(&BinaryPoint::Infinity)
                    .unwrap();
                assert!(infinity.identity_holds);

                for u in 0..8 {
                    for v in 0..8 {
                        let point = BinaryCoverPoint::Affine {
                            u: elt(u),
                            v: elt(v),
                        };
                        if !c.source_is_on_curve(&point) {
                            continue;
                        }
                        source_points += 1;
                        let image = c.pushforward_rational_point_divisor(&point).unwrap();
                        assert!(c.target_is_on_curve(&image));
                        if u != 0 {
                            let pullback = c.pullback_point_class(&image).unwrap();
                            assert!(pullback.u.eval(&elt(u), &c.irreducible).is_zero());
                            assert_eq!(pullback.v.eval(&elt(u), &c.irreducible), elt(v));
                        }
                    }
                }
            }
        }
        assert!(target_points > 0);
        assert!(source_points > 0);
        assert!(generic_pullbacks > 0);
    }

    #[test]
    fn subgroup_status_and_arbitrary_pushforward_boundary() {
        let c = cover(0, 1);
        assert_eq!(
            c.subgroup_transfer_status(&BigUint::from(7u32)).unwrap(),
            SubgroupTransferStatus::InjectiveIfExactSubgroupOrder {
                order: BigUint::from(7u32),
                degree_inverse_mod_order: BigUint::from(5u32),
            }
        );
        assert_eq!(
            c.subgroup_transfer_status(&BigUint::from(21u32)).unwrap(),
            SubgroupTransferStatus::KernelIntersectionUnresolved {
                order: BigUint::from(21u32),
                common_factor: 3,
            }
        );
        assert_eq!(
            c.subgroup_transfer_status(&BigUint::zero()).unwrap_err(),
            CoverTransferError::InvalidSubgroupOrder
        );
        let id = MumfordDivisor::identity(3);
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&id).unwrap_err(),
            CoverTransferError::Unsupported(
                UnsupportedCoverOperation::ArbitraryJacobianPushforward
            )
        );

        // Exact EC1/ICV1 catalog fixture:
        // EC1N7Ce0hb6d297a2ca08 on z^7+z+1, subgroup order 29,
        // generator (0x26,0x5).  Because 29 is prime, [29]G=O and G!=O
        // certify the order.  The checked pullback has the same order, and
        // multiplying its pushforward [3]G by 3^{-1}=10 recovers G.
        let catalog = OrdinaryBinaryCover::new(
            7,
            IrreduciblePoly {
                degree: 7,
                low_terms: vec![0, 1],
            },
            F2mElement::zero(7),
            F2mElement::one(7),
        )
        .unwrap();
        let generator = BinaryPoint::Affine {
            x: F2mElement::from_biguint(&BigUint::from(0x26u32), 7),
            y: F2mElement::from_biguint(&BigUint::from(0x05u32), 7),
        };
        assert!(catalog.target_is_on_curve(&generator));
        assert_ne!(generator, BinaryPoint::Infinity);
        assert_eq!(
            catalog.target_scalar_mul_small(&generator, 29).unwrap(),
            BinaryPoint::Infinity
        );
        assert_eq!(
            catalog
                .subgroup_transfer_status(&BigUint::from(29u32))
                .unwrap(),
            SubgroupTransferStatus::InjectiveIfExactSubgroupOrder {
                order: BigUint::from(29u32),
                degree_inverse_mod_order: BigUint::from(10u32),
            }
        );
        let image = catalog.pullback_point_class(&generator).unwrap();
        assert!(!image.is_identity());
        assert!(image
            .try_scalar_mul(&BigUint::from(29u32), catalog.source_curve())
            .unwrap()
            .is_identity());
        let composition = catalog
            .pull_push_composition_certificate(&generator)
            .unwrap();
        assert_eq!(
            catalog
                .target_scalar_mul_small(&composition.formal_pushforward, 10)
                .unwrap(),
            generator
        );
    }

    #[test]
    fn fixed_f8_arbitrary_pushforward_is_complete_and_homomorphic() {
        let c = cover(1, 3);
        let divisors = all_fixed_f8_reduced_divisors(c.source_curve());
        let divisor_keys: HashSet<_> = divisors.iter().map(fixed_f8_divisor_key).collect();
        let mut degree_counts = [0usize; 4];
        for divisor in &divisors {
            divisor.validate(c.source_curve()).unwrap();
            degree_counts[divisor.u.degree().unwrap()] += 1;
        }
        assert_eq!(degree_counts, [1, 5, 56, 424]);
        assert_eq!(divisors.len(), 486);
        assert_eq!(divisor_keys.len(), 486);

        let identity = MumfordDivisor::identity(3);
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&identity).unwrap(),
            BinaryPoint::Infinity
        );
        let table = FIXED_F8_NORM_TABLE.get().unwrap().as_ref().unwrap();
        assert_eq!(table.images.len(), 486);
        assert_eq!(
            table.images.keys().cloned().collect::<HashSet<_>>(),
            divisor_keys
        );

        let q = target_point(4, 4);
        let target_group = c.fixed_f8_target_multiples(&q, 6).unwrap();
        assert_eq!(target_group.len(), 6);
        assert_eq!(
            c.target_scalar_mul_small(&q, 6).unwrap(),
            BinaryPoint::Infinity
        );
        let mut fiber_sizes = [0usize; 6];
        let mut images = Vec::with_capacity(divisors.len());
        for divisor in &divisors {
            let image = c.pushforward_arbitrary_jacobian_class(divisor).unwrap();
            assert!(c.target_is_on_curve(&image));
            let image_index = target_group
                .iter()
                .position(|candidate| candidate == &image)
                .expect("fixed-fixture norm image lies in the six-element target");
            fiber_sizes[image_index] += 1;

            let inverse = divisor.try_neg(c.source_curve()).unwrap();
            let inverse_image = c.pushforward_arbitrary_jacobian_class(&inverse).unwrap();
            assert_eq!(
                c.target_add(&image, &inverse_image).unwrap(),
                BinaryPoint::Infinity
            );
            images.push(image);
        }
        assert_eq!(fiber_sizes, [81; 6]);

        for (left_index, left) in divisors.iter().enumerate() {
            for (right_index, right) in divisors.iter().enumerate() {
                let sum = left.try_add(right, c.source_curve()).unwrap();
                let actual = c.pushforward_arbitrary_jacobian_class(&sum).unwrap();
                let expected = c
                    .target_add(&images[left_index], &images[right_index])
                    .unwrap();
                assert_eq!(
                    actual, expected,
                    "norm homomorphism failed for ordered pair ({left_index},{right_index})"
                );
            }
        }

        let mut rational_source_points = vec![BinaryCoverPoint::Infinity];
        for u in 0..8 {
            for v in 0..8 {
                let point = BinaryCoverPoint::Affine {
                    u: elt(u),
                    v: elt(v),
                };
                if c.source_is_on_curve(&point) {
                    rational_source_points.push(point);
                }
            }
        }
        assert_eq!(rational_source_points.len(), 6);
        for point in rational_source_points {
            let divisor = match &point {
                BinaryCoverPoint::Infinity => MumfordDivisor::identity(3),
                BinaryCoverPoint::Affine { u, v } => {
                    MumfordDivisor::try_from_point(c.source_curve(), u, v).unwrap()
                }
            };
            assert_eq!(
                c.pushforward_arbitrary_jacobian_class(&divisor).unwrap(),
                c.pushforward_rational_point_divisor(&point).unwrap()
            );
        }

        let mut target_points = vec![BinaryPoint::Infinity];
        for x in 0..8 {
            for y in 0..8 {
                let point = target_point(x, y);
                if c.target_is_on_curve(&point) {
                    target_points.push(point);
                }
            }
        }
        assert_eq!(target_points.len(), 6);
        for point in target_points {
            let pullback = c.pullback_point_class(&point).unwrap();
            assert_eq!(
                c.pushforward_arbitrary_jacobian_class(&pullback).unwrap(),
                c.target_scalar_mul_small(&point, 3).unwrap()
            );
        }

        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&checked_divisor(
                c.source_curve(),
                &[0, 1],
                &[],
            ))
            .unwrap(),
            target_point(0, 7)
        );
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&checked_divisor(
                c.source_curve(),
                &[4, 0, 1],
                &[2, 5],
            ))
            .unwrap(),
            target_point(3, 2)
        );
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&checked_divisor(
                c.source_curve(),
                &[1, 1, 1],
                &[3, 2],
            ))
            .unwrap(),
            BinaryPoint::Infinity
        );
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&checked_divisor(
                c.source_curve(),
                &[6, 1, 0, 1],
                &[6, 7, 6],
            ))
            .unwrap(),
            target_point(4, 0)
        );
    }

    #[test]
    fn arbitrary_pushforward_rejects_malformed_and_unsupported_inputs() {
        let c = cover(1, 3);
        let identity = MumfordDivisor::identity(3);
        let zero_u = MumfordDivisor {
            u: F2mPoly::zero(3),
            v: F2mPoly::zero(3),
        };
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&zero_u),
            Err(CoverTransferError::Hyperelliptic(
                HyperellipticError::ZeroDivisorPolynomial
            ))
        );
        let non_monic = MumfordDivisor {
            u: F2mPoly::constant(elt(2)),
            v: F2mPoly::zero(3),
        };
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&non_monic),
            Err(CoverTransferError::Hyperelliptic(
                HyperellipticError::NonMonicDivisorPolynomial
            ))
        );
        let wrong_width = MumfordDivisor::identity(2);
        assert!(matches!(
            c.pushforward_arbitrary_jacobian_class(&wrong_width),
            Err(CoverTransferError::Hyperelliptic(
                HyperellipticError::PolynomialFieldMismatch { .. }
            ))
        ));
        let wrong_equation = MumfordDivisor {
            u: F2mPoly::from_coeffs(vec![elt(1), elt(1)], 3),
            v: F2mPoly::zero(3),
        };
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&wrong_equation),
            Err(CoverTransferError::Hyperelliptic(
                HyperellipticError::DivisorEquationNotDivisible
            ))
        );

        for a in 0..8 {
            for b in 1..8 {
                if (a, b) == (1, 3) {
                    continue;
                }
                assert_eq!(
                    cover(a, b).pushforward_arbitrary_jacobian_class(&identity),
                    Err(CoverTransferError::Unsupported(
                        UnsupportedCoverOperation::ArbitraryJacobianPushforward
                    )),
                    "model a={a}, b={b} must remain outside the bounded fixture"
                );
            }
        }
        let other_basis = OrdinaryBinaryCover::new(
            3,
            IrreduciblePoly {
                degree: 3,
                low_terms: vec![0, 2],
            },
            elt(1),
            elt(3),
        )
        .unwrap();
        assert_eq!(
            other_basis.pushforward_arbitrary_jacobian_class(&identity),
            Err(CoverTransferError::Unsupported(
                UnsupportedCoverOperation::ArbitraryJacobianPushforward
            ))
        );
    }

    #[test]
    fn arbitrary_pushforward_rejects_oversized_polynomials_before_validation() {
        let c = cover(1, 3);
        let oversized_u = MumfordDivisor {
            u: F2mPoly {
                m: 2,
                coeffs: vec![F2mElement::one(2); 5],
            },
            v: F2mPoly::zero(3),
        };
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&oversized_u),
            Err(CoverTransferError::DivisorPolynomialTooLong {
                polynomial: "u",
                maximum_coefficients: 4,
                actual_coefficients: 5,
            })
        );

        let oversized_v = MumfordDivisor {
            u: F2mPoly::one(3),
            v: F2mPoly {
                m: 2,
                coeffs: vec![F2mElement::one(2); 4],
            },
        };
        assert_eq!(
            c.pushforward_arbitrary_jacobian_class(&oversized_v),
            Err(CoverTransferError::DivisorPolynomialTooLong {
                polynomial: "v",
                maximum_coefficients: 3,
                actual_coefficients: 4,
            })
        );
    }
}
