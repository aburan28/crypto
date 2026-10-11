//! Research-only helpers for the Eigen-Séta design note.
//!
//! This module records and tests the algebra that is needed by the proposed
//! construction.  It is **not** a public-key encryption implementation:
//! key generation, serialization, constant-time arithmetic, random-curve
//! generation, and IND-CCA security are intentionally not provided.
//!
//! The module has three deliberately small responsibilities:
//!
//! * validate the inequalities used by the toy parameter discussion;
//! * model the diagonal mask / eigenbasis commutation lemma; and
//! * provide a typed wrapper for the existing single `(N,N)` theta chain.
//!
//! Keeping these pieces separate is important: a passing Kani split test is a
//! correctness result for an algebraic subroutine, not evidence that the
//! Eigen-Séta encryption proposal is secure.

use crate::curve::{Curve, Pt};
use crate::field::Field;
use crate::theta::{chain, ChainResult};

/// Toy parameter set for the design note.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Parameters {
    /// Security parameter used only as a bookkeeping label.
    pub lambda: u32,
    /// The sender isogeny degree is represented as `3^degree_exponent`.
    pub degree_exponent: u32,
    /// The torsion order is represented as `2^torsion_exponent`.
    pub torsion_exponent: u32,
    /// Absolute upper bound on the toy field prime.
    pub prime: u64,
    /// Absolute value of the trace-zero endomorphism norm in the toy model.
    pub delta: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ParameterError {
    ZeroParameter,
    Overflow,
    TorsionDoesNotDominate,
    PrimeDoesNotDominate,
}

impl Parameters {
    /// Returns `3^degree_exponent`, when it fits in `u64`.
    pub fn sender_degree(self) -> Result<u64, ParameterError> {
        if self.degree_exponent == 0 {
            return Err(ParameterError::ZeroParameter);
        }
        checked_pow(3, self.degree_exponent)
    }

    /// Returns `2^torsion_exponent`, when it fits in `u64`.
    pub fn torsion_order(self) -> Result<u64, ParameterError> {
        if self.torsion_exponent == 0 {
            return Err(ParameterError::ZeroParameter);
        }
        checked_pow(2, self.torsion_exponent)
    }

    /// Checks the toy analogue of `N > d * sqrt(delta)` and `p > N*d`.
    ///
    /// The comparison is squared to avoid floating point arithmetic:
    /// `N² > d²*delta`, followed by `p > N*d`.
    pub fn validate(self) -> Result<(), ParameterError> {
        if self.lambda == 0 || self.prime == 0 || self.delta == 0 {
            return Err(ParameterError::ZeroParameter);
        }
        let d = self.sender_degree()? as u128;
        let n = self.torsion_order()? as u128;
        let n2 = n.checked_mul(n).ok_or(ParameterError::Overflow)?;
        let d2_delta = d
            .checked_mul(d)
            .and_then(|value| value.checked_mul(self.delta as u128))
            .ok_or(ParameterError::Overflow)?;
        if n2 <= d2_delta {
            return Err(ParameterError::TorsionDoesNotDominate);
        }
        let nd = n.checked_mul(d).ok_or(ParameterError::Overflow)?;
        if self.prime as u128 <= nd {
            return Err(ParameterError::PrimeDoesNotDominate);
        }
        Ok(())
    }
}

fn checked_pow(base: u64, exponent: u32) -> Result<u64, ParameterError> {
    let mut out = 1u64;
    for _ in 0..exponent {
        out = out.checked_mul(base).ok_or(ParameterError::Overflow)?;
    }
    Ok(out)
}

/// A diagonal mask in a published eigenbasis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DiagonalMask<E> {
    pub first: E,
    pub second: E,
}

/// A diagonal action of the secret trace-zero endomorphism on the eigenbasis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EigenAction<E> {
    pub first: E,
    pub second: E,
}

impl<E: Copy + PartialEq> DiagonalMask<E> {
    /// The commutation condition used by the lollipop decoding lemma.
    pub fn commutes_with(self, _action: EigenAction<E>) -> bool {
        // Both matrices are diagonal.  The explicit method keeps the
        // centralizer claim visible at call sites and leaves room for a
        // future non-diagonal action type without changing the transcript API.
        true
    }

    /// Applies the masked diagonal action to a two-coordinate torsion image.
    pub fn apply(self, action: EigenAction<E>, point: [E; 2], mul: impl Fn(E, E) -> E) -> [E; 2] {
        [
            mul(self.first, mul(action.first, point[0])),
            mul(self.second, mul(action.second, point[1])),
        ]
    }
}

/// A non-secret request for one `(N,N)` Kani chain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KaniRequest<E> {
    /// The two lifted kernel generators, each living on a product of curves.
    pub generators: [(Pt<E>, Pt<E>); 2],
    /// Number of `(2,2)` steps after the initial gluing.
    pub steps: u32,
}

/// Runs one `(N,N)` theta chain.  This does not compute the second half of
/// the `(N²,N²)` isogeny required by the Eigen-Séta design note.  Callers
/// must construct and validate all torsion points.
pub fn evaluate_kani<F: Field>(
    field: &F,
    left: &Curve<F::E>,
    right: &Curve<F::E>,
    request: KaniRequest<F::E>,
    extra: &[(Pt<F::E>, Pt<F::E>)],
) -> Option<ChainResult<F::E>> {
    if request.steps == 0 {
        return None;
    }
    chain(field, left, right, request.generators, request.steps, extra)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toy_parameter_checks_are_fail_closed() {
        let good = Parameters {
            lambda: 4,
            degree_exponent: 2,
            torsion_exponent: 8,
            prime: 10_000,
            delta: 7,
        };
        assert_eq!(good.validate(), Ok(()));

        let bad = Parameters { prime: 100, ..good };
        assert_eq!(bad.validate(), Err(ParameterError::PrimeDoesNotDominate));
    }

    #[test]
    fn diagonal_masks_apply_coordinatewise() {
        let mask = DiagonalMask {
            first: 3u64,
            second: 5u64,
        };
        let action = EigenAction {
            first: 7u64,
            second: 11u64,
        };
        assert!(mask.commutes_with(action));
        assert_eq!(mask.apply(action, [2, 4], |a, b| a * b), [42, 220]);
    }
}
