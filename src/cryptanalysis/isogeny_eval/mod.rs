//! **Isogeny evaluation over prime fields `F_p`, `5 ≤ p < 2^63`**, with
//! optional operation counting.
//!
//! Moved out of the round-4 isogeny-fibers examples
//! (`examples/r4_16_3_vertical_margin.rs`,
//! `examples/r4_16_5_class_multitarget_rho.rs`, which stay frozen as the
//! evidence for `RESEARCH_ISOGENY_FIBERS_ECDLP.md` §17.3 and §17.5).
//!
//! * [`counted_fp`]: `F_p` on `u64` words, generic over an
//!   [`counted_fp::OpCounter`].  [`counted_fp::PrimeField::new`] counts
//!   nothing (a zero-sized counter, free); [`counted_fp::PrimeField::counted`]
//!   counts multiplications `M` and inversions `I`.
//! * [`poly`]: Karatsuba multiplication, Newton reciprocal, remainder by a
//!   monic divisor, balanced products, product and remainder trees.
//! * [`montgomery`]: x-only Montgomery `xDBL`, `xADD`, ladder (generic over
//!   [`montgomery::FieldArith`]), batch normalisation, `j`.
//! * [`velu`]: `O(ℓ)` Vélu for odd prime `ℓ` (Costello–Hisil images,
//!   Meyer–Reith codomain).
//! * [`sqrt_velu`]: √élu (Bernstein–De Feo–Leroux–Smith 2020).
//! * [`small_degree`]: `ℓ ∈ {2, 3}` on short Weierstrass curves with `(x, y)`
//!   point evaluation, `F_p`-isomorphisms and a certified dual.
//!
//! Scope and limits: prime fields only (no `F_{p^k}` type; the x-only
//! formulas are generic so one can be added), `p < 2^63`, kernels generated
//! by an `F_p`-rational x-coordinate.  Polynomial multiplication is
//! Karatsuba; NTT/FFT multiplication (needed for the asymptotic `Õ(√ℓ)` of
//! √élu) is a follow-up.  The existing library pieces this complements:
//! `crate::cryptanalysis::isogeny_walk::kernel` (odd-degree codomains from
//! a kernel polynomial over a 256-bit field, used by the tests here as an
//! independent check), `crate::isogeny::velu` (brute-force toy Vélu),
//! `crate::pqc::fast::isogeny` (2-isogeny chains over `F_{p²}`), and
//! `crate::cryptanalysis::binary_velu` (characteristic 2).

pub mod counted_fp;
pub mod montgomery;
pub mod poly;
pub mod small_degree;
pub mod sqrt_velu;
pub mod velu;

#[cfg(test)]
mod tests;

pub use counted_fp::{NoCount, OpCounter, OpCounts, OpTally, PrimeField};
pub use montgomery::{MontgomeryCurve, XZ};
pub use small_degree::{ShortWeierstrass, SmallIsogeny, WPoint};
pub use sqrt_velu::{sqrt_velu, SqrtVeluShape};
pub use velu::{velu, IsogenyEvaluation};

/// Why an isogeny could not be built or evaluated.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsogenyEvalError {
    /// The degree is not supported here (not an odd prime for
    /// [`velu()`]/[`sqrt_velu()`], not 2 or 3 for [`SmallIsogeny`]).
    InvalidDegree(u64),
    /// A curve (given or computed) is singular.
    SingularCurve,
    /// The given kernel x-coordinate does not generate a kernel of the
    /// stated degree.
    NotAKernelPoint,
    /// No point of large enough order was found to certify a dual.
    NoCertificationPoint,
    /// No candidate dual passed the certificate.
    NoDual,
}

impl std::fmt::Display for IsogenyEvalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IsogenyEvalError::InvalidDegree(l) => write!(f, "unsupported isogeny degree {l}"),
            IsogenyEvalError::SingularCurve => write!(f, "singular curve"),
            IsogenyEvalError::NotAKernelPoint => {
                write!(f, "x-coordinate does not generate a kernel of that degree")
            }
            IsogenyEvalError::NoCertificationPoint => {
                write!(f, "no point of order > 4ℓ² found to certify the dual")
            }
            IsogenyEvalError::NoDual => write!(f, "no certified dual isogeny found"),
        }
    }
}

impl std::error::Error for IsogenyEvalError {}
