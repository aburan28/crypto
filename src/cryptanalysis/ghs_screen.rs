//! Native structural screener for classical binary GHS Weil descent.
//!
//! Given an ordinary binary curve
//! `E: y^2 + x*y = x^3 + a*x^2 + b` over a polynomial-basis field
//! `F_{2^N}`, this module validates the field and curve input, enumerates every
//! non-trivial factorisation `N = n*l`, and reports the exact GHS magic number,
//! genus, and Artin--Schreier cover degree for the tower
//! `F_2 <= F_{2^l} <= F_{2^N}`.
//!
//! This is a structural screen.  A small genus is a candidate for further
//! analysis, not an end-to-end attack result: subgroup transport, Jacobian
//! order, relation collection, linear algebra, and comparison with a matched
//! Pollard-rho reference are deliberately outside this module.

use crate::{
    binary_ecc::{F2mElement, IrreduciblePoly},
    cryptanalysis::ec_trapdoor::audit_curve,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::{collections::HashSet, error::Error, fmt};

/// Largest accepted field degree.  This matches the curve-cover catalog's
/// native checker and prevents hostile inputs from requesting unbounded field
/// allocations.
pub const MAX_FIELD_DEGREE: u32 = 4096;

/// An ordinary binary curve and its exact polynomial-basis field.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GhsCurveInput {
    /// Absolute degree `N` in `F_{2^N}`.
    pub absolute_degree: u32,
    /// Monic binary polynomial, including its `z^N` bit.
    pub modulus: BigUint,
    /// Curve coefficient `a`, encoded as a polynomial bitset.
    pub a: BigUint,
    /// Nonzero curve coefficient `b`, encoded as a polynomial bitset.
    pub b: BigUint,
}

/// One exact GHS cover tower attached to a field factorisation `N = n*l`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GhsCoverTower {
    /// `[F_{2^N}:F_{2^l}]`.
    pub relative_degree: u32,
    /// Base-field degree `l`.
    pub base_degree: u32,
    /// Hess magic number for this factorisation.
    pub magic_number: u32,
    /// Exact GHS descent genus (never truncated to a machine word).
    pub genus: BigUint,
    /// Whether the existing GHS audit selected its type-I branch.
    pub type_i: bool,
    /// Degree `2^magic_number` of the Artin--Schreier compositum cover.
    pub cover_degree: BigUint,
    /// Whether this row satisfies the caller's structural genus bound.
    pub within_genus_bound: bool,
}

impl GhsCoverTower {
    /// Human-readable constant/base/top field chain.
    pub fn field_tower(&self, absolute_degree: u32) -> String {
        format!(
            "F_2 <= F_(2^{}) <= F_(2^{absolute_degree})",
            self.base_degree
        )
    }

    /// Human-readable curve side of the GHS construction.
    pub fn cover_tower(&self, absolute_degree: u32) -> String {
        format!("C/F_(2^{}) -> E/F_(2^{absolute_degree})", self.base_degree)
    }
}

/// Complete, genus-sorted structural screen.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GhsScreenReport {
    pub absolute_degree: u32,
    pub genus_bound: BigUint,
    pub rows: Vec<GhsCoverTower>,
}

impl GhsScreenReport {
    /// Lowest-genus factorisation, if the absolute degree has one.
    pub fn best(&self) -> Option<&GhsCoverTower> {
        self.rows.first()
    }

    /// Whether at least one factorisation is inside the requested genus bound.
    pub fn has_candidate(&self) -> bool {
        self.rows.iter().any(|row| row.within_genus_bound)
    }
}

/// Invalid field or curve supplied to [`screen_curve`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GhsScreenError(String);

impl GhsScreenError {
    fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}

impl fmt::Display for GhsScreenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl Error for GhsScreenError {}

/// Validate and screen an ordinary binary curve across every factorisation of
/// its absolute field degree.
pub fn screen_curve(
    input: &GhsCurveInput,
    genus_bound: &BigUint,
) -> Result<GhsScreenReport, GhsScreenError> {
    let irr = validate_input(input)?;
    let a = F2mElement::from_biguint(&input.a, input.absolute_degree);
    let b = F2mElement::from_biguint(&input.b, input.absolute_degree);
    let rows = audit_curve(input.absolute_degree, &irr, &a, &b)
        .into_iter()
        .map(|row| {
            let cover_degree = BigUint::one() << row.magic_m;
            let within_genus_bound = &row.genus <= genus_bound;
            GhsCoverTower {
                relative_degree: row.n,
                base_degree: row.l,
                magic_number: row.magic_m,
                genus: row.genus,
                type_i: row.type_i,
                cover_degree,
                within_genus_bound,
            }
        })
        .collect();
    Ok(GhsScreenReport {
        absolute_degree: input.absolute_degree,
        genus_bound: genus_bound.clone(),
        rows,
    })
}

fn validate_input(input: &GhsCurveInput) -> Result<IrreduciblePoly, GhsScreenError> {
    let degree = input.absolute_degree;
    if !(2..=MAX_FIELD_DEGREE).contains(&degree) {
        return Err(GhsScreenError::new(format!(
            "absolute degree must be in 2..={MAX_FIELD_DEGREE}"
        )));
    }
    if input.modulus.bits() != u64::from(degree) + 1 || !input.modulus.bit(0) {
        return Err(GhsScreenError::new(
            "modulus must be monic of the declared degree with nonzero constant term",
        ));
    }
    if input.a.bits() > u64::from(degree) || input.b.bits() > u64::from(degree) {
        return Err(GhsScreenError::new(
            "curve coefficients must be canonical field elements",
        ));
    }
    if input.b.is_zero() {
        return Err(GhsScreenError::new(
            "ordinary binary model requires nonzero b",
        ));
    }
    let irr = IrreduciblePoly {
        degree,
        low_terms: (0..degree)
            .filter(|exponent| input.modulus.bit(u64::from(*exponent)))
            .collect(),
    };
    if !is_irreducible(&input.modulus, &irr) {
        return Err(GhsScreenError::new("binary modulus is reducible"));
    }
    Ok(irr)
}

fn polynomial_remainder(mut dividend: BigUint, divisor: &BigUint) -> BigUint {
    while !dividend.is_zero() && dividend.bits() >= divisor.bits() {
        dividend ^= divisor << ((dividend.bits() - divisor.bits()) as usize);
    }
    dividend
}

fn polynomial_gcd(mut left: BigUint, mut right: BigUint) -> BigUint {
    while !right.is_zero() {
        let remainder = polynomial_remainder(left, &right);
        left = right;
        right = remainder;
    }
    left
}

fn prime_division_checkpoints(degree: u32) -> HashSet<u32> {
    let mut quotient = degree;
    let mut prime = 2u32;
    let mut checkpoints = HashSet::new();
    while prime * prime <= quotient {
        if quotient.is_multiple_of(prime) {
            checkpoints.insert(degree / prime);
            while quotient.is_multiple_of(prime) {
                quotient /= prime;
            }
        }
        prime += 1;
    }
    if quotient > 1 {
        checkpoints.insert(degree / quotient);
    }
    checkpoints
}

/// Rabin's exact irreducibility criterion over `F_2`.
fn is_irreducible(modulus: &BigUint, irr: &IrreduciblePoly) -> bool {
    let z = polynomial_remainder(BigUint::from(2u32), modulus);
    let checkpoints = prime_division_checkpoints(irr.degree);
    let mut power = F2mElement::from_biguint(&z, irr.degree);
    for exponent in 1..=irr.degree {
        power = power.square(irr);
        if checkpoints.contains(&exponent)
            && polynomial_gcd(power.to_biguint() ^ &z, modulus.clone()) != BigUint::one()
        {
            return false;
        }
    }
    power.to_biguint() == z
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aes_curve(b: u32) -> GhsCurveInput {
        GhsCurveInput {
            absolute_degree: 8,
            modulus: BigUint::from(0x11bu32),
            a: BigUint::zero(),
            b: BigUint::from(b),
        }
    }

    #[test]
    fn enumerates_every_factorisation_and_cover_tower() {
        let report = screen_curve(&aes_curve(1), &BigUint::from(4u32)).unwrap();
        let mut factors: Vec<_> = report
            .rows
            .iter()
            .map(|row| (row.relative_degree, row.base_degree))
            .collect();
        factors.sort_unstable();
        assert_eq!(factors, vec![(2, 4), (4, 2), (8, 1)]);
        assert!(report.has_candidate());
        for row in &report.rows {
            assert_eq!(row.cover_degree, BigUint::from(2u32));
            assert_eq!(
                row.field_tower(8),
                format!("F_2 <= F_(2^{}) <= F_(2^8)", row.base_degree)
            );
            assert_eq!(
                row.cover_tower(8),
                format!("C/F_(2^{}) -> E/F_(2^8)", row.base_degree)
            );
        }
    }

    #[test]
    fn folds_a_contribution_into_each_rows_genus() {
        let mut input = aes_curve(1);
        input.a = BigUint::from(2u32);
        let report = screen_curve(&input, &BigUint::from(4u32)).unwrap();
        let rows: Vec<_> = report
            .rows
            .iter()
            .map(|row| {
                (
                    row.relative_degree,
                    row.magic_number,
                    row.genus.clone(),
                    row.type_i,
                )
            })
            .collect();
        assert_eq!(
            rows,
            vec![
                (2, 2, BigUint::from(2u32), false),
                (4, 3, BigUint::from(4u32), false),
                (8, 6, BigUint::from(32u32), false),
            ]
        );
    }

    #[test]
    fn rejects_reducible_fields_and_singular_models() {
        let mut reducible = aes_curve(1);
        reducible.modulus = BigUint::from(0x101u32); // z^8 + 1 = (z + 1)^8.
        assert_eq!(
            screen_curve(&reducible, &BigUint::from(4u32))
                .unwrap_err()
                .to_string(),
            "binary modulus is reducible"
        );
        assert_eq!(
            screen_curve(&aes_curve(0), &BigUint::from(4u32))
                .unwrap_err()
                .to_string(),
            "ordinary binary model requires nonzero b"
        );
    }
}
